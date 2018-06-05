clear; clc; close all;
addpath('../', '../urdf/')
save_dir = ['./data_28-May-2018' , '/'];
trials = dir([save_dir, '*TROT*.mat']);

global lp_b

lp_b = [0, 7.58, -11.350;
    0, 7.58, -11.350;
    0, -7.58, -11.350;
    0, -7.58, -11.350];

DEBUG = 0;

%% Load Rigid Body Model

urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', 'HAMR_scaledV2.urdf');
dt = 0.4;

% general rbm options
options.floating = true;
options.collision = true;
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.z_inactive_guess_tol = 0.1;
options.use_bullet = false;
options.terrain = RigidBodyFlatTerrain();

% hamr options
options.stiffness_mult = [1, 1, 1, 1]';     %(swing act, lift act, swing flex, lift flex)

% Build robot + visualizer
hamr = HamrTSRBM(HamrRBM(urdf, options), dt, options);
nq = hamr.getNumPositions();
nv = hamr.getNumVelocities();
nu = hamr.getNumInputs();
qa = hamr.getActuatedJoints();

v = hamr.constructVisualizer();

%%
Nf = numel(trials); % size(cl_sols, 1);
data0 = load([trials(1).folder, '/', trials(1).name]);
Ndc_swing = size(data0.cl_sols, 1);
Np_lift = size(data0.cl_sols, 2);
DC_SWING = zeros(Ndc_swing, Np_lift);
P_LIFT = zeros(Ndc_swing, Np_lift);
freq = zeros(Nf, 1);

avg_speed = zeros(Ndc_swing, Np_lift, Nf);
avg_slip = zeros(Ndc_swing, Np_lift, Nf);
% rms_err = zeros(Np_lift, Ndc_swing, Nf);

% fixed parameters
VB = 225;

for i = 1:Nf % frequency
    
    datai = load([trials(i).folder, '/',  trials(i).name]);
    disp(trials(i).name)
    cl_sols = datai.cl_sols;
    
    for j = 1:Ndc_swing % swing duty cycle (rows)
        for m = 1:Np_lift % lift push (columns)
            
            % unpack parameters
            sol_struct = cl_sols{j, m};
            freq(i) = sol_struct.params.FREQ;
            DC_SWING(j,m) = sol_struct.params.DC_SWING;
            P_LIFT(j,m) = sol_struct.params.DC_LIFT;
            RAMPCYC = sol_struct.params.RAMPCYC;
            TOTCYC = sol_struct.params.NCYC;
            DATACYC = TOTCYC - RAMPCYC;
            
            
            fprintf('(%d, %d) is %d swing and %d lift \r', j, m, ...
                DC_SWING(j,m),  P_LIFT(j,m))
            
            
            % compute average speed
            tt_sol = sol_struct.tt_sol;
            dt_sol = mean(diff(tt_sol));
            xx_sol = sol_struct.xx_sol;
            ind0 = find(tt_sol > RAMPCYC/freq(i), 1, 'first');
            x = xx_sol(1,ind0:end) - xx_sol(1,ind0);
            y = xx_sol(2,ind0:end) - xx_sol(2,ind0);
            avg_speed(j,m,i) = sum(sqrt(diff(x).^2 + diff(y).^2))/(tt_sol(end) - tt_sol(ind0));
            %         avg_speed(j,m,i) = mean(xx_sol(39, ind0:end)); %(xx_sol(1, end) - xx_sol(1, ind0))/(tt_sol(end) - tt_sol(ind0));
            %         avg_speed(j,m,i) = (xx_sol(1, end) - xx_sol(1, ind0))/(tt_sol(end) - tt_sol(ind0));
            
            if DEBUG
                xtraj = PPTrajectory(foh(tt_sol/1e3, xx_sol(1:nq,:)));
                xtraj = xtraj.setOutputFrame(v.getInputFrame());
                v.playback(xtraj, struct('slider', true));
            end
            %
            figure(i); hold on;
            subplot(2,1,1); hold on;
            title([num2str(1e3*freq(i)) ' Hz'])
            plot(tt_sol, xx_sol(nq+1, :));
            subplot(2,1,2); hold on;
            plot(tt_sol, xx_sol(nq+2, :));
            %         legend_str{Ndc_swing*(j-1) + m} = [num2str(DC_SWING(m,j)) '-' num2str(P_LIFT(m,j))];
            
            % compute slip
            pfCL = sol_struct.pfCL;
            vfCL = sol_struct.vfCL;
            
            zmin = min(pfCL(:,3,:));
            vfX = vfCL(:,1,:);
            vfY = vfCL(:,2,:);
            xfTot = trapz(tt_sol, abs(vfX));
            
            xfB = zeros(nu/2, 1);
            yF = zeros(nu/2, 1);
            slipf = 0;
            for f =1:(nu/2)
                if sum(tt_sol(vfX(:,1,f)<0)) == 0
                    xfB(f) = 0;
                else
                    xfB(f) = sum(abs(vfX(vfX(:,1,f)<0, 1, f))*dt_sol);
                    z0mask = pfCL(1:end-1, 3, f)<zmin(1,1,f);
                    yF(f) =  sum(abs(vfY(z0mask, 1, f))*dt_sol);
                end
                slipf = slipf + (1 - (xfB(f)+yF(f))/xfTot(1,1,f));
            end
            avg_slip(j,m,i) = slipf/(nu/2);
        end
    end
    %     legend(legend_str)
end

%%
% figure(10); hold on;
% plot(freq, max(max(avg_speed,[], 3), [], 2), '*')
% plot(freq, mean(mean(avg_speed, 3), 2), '+')
%
% plot(freq, 10*freq)
% plot(freq, 5*freq)

for i = 1:Nf
    figure(Nf+i); clf;
    subplot(2,1,1); hold on;
    title([num2str(1e3*freq(i)) ' Hz'])
    contourf(DC_SWING, P_LIFT, avg_speed(:,:,i)/(freq(i)))
    xlabel('swing dc')
    ylabel('lift dc')
    colorbar
    subplot(2,1,2); hold on;
    title([num2str(1e3*freq(i)) ' Hz'])
    contourf(DC_SWING, P_LIFT, avg_slip(:,:,i))
    colorbar
    
end
%
tilefigs;
