clear; clc; close all;
warning('off', 'MATLAB:nargchk:deprecated');
addpath('../', '../urdf/')

save_dir = ['./data_' datestr(date), '/'];
if exist(save_dir, 'dir') == 0
    mkdir(save_dir);
end

% Build robot + visualizer
urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', 'HAMR_scaledV2.urdf');
options.floating = false;
options.collision = false;
hamrID = HamrRBMID(urdf, options);
v = hamrID.constructVisualizer();

% build actuators
dp.Vb = 225;
dp.Vg = 0;

nact = 8;
actuators = HamrActuators(nact, {'FLsact', 'FLlact', 'RLsact', 'RLlact', ...
    'FRsact', 'FRlact', 'RRsact', 'RRlact'}, [], dp);

% dimensions
nq = hamrID.getNumPositions();
nv = hamrID.getNumVelocities();
nu = hamrID.getNumInputs();
nl = hamrID.NL;
qa = hamrID.getActuatedJoints();
nqa = numel(qa);

% build desired actuator trajectory
gait = 'TROT';
freq = linspace(10, 50, 5)*1e-3;        % frequency
DC = linspace(50, 80, 4);                  % duty cycle for swing
DL = linspace(20, 75, 5);                  % percent "push" into the ground
LIFTAMP = 0.15;                                 % lift actuator motion (mm)
SWINGAMP = 0.175;                               % swing actuator motion (mm)
NSAMP = 100;                                    % number of samples
TYPE = 1; 

Nf = numel(freq);
Ndc = numel(DC);
Ndl = numel(DL);

for i = 1:Nf
    ffinputs = cell(Ndc, Ndl);
    
    for j = 1:Ndc
        for k = 1:Ndl
            
            [trajd, brkVals] = GenerateHueristicActTraj(gait, freq(i), DC(j), DL(k), NSAMP, TYPE);
            
            trajd(:, 2:2:end) = SWINGAMP*trajd(:, 2:2:end);
            trajd(:, 3:2:end) = LIFTAMP*trajd(:, 3:2:end);
%             trajd = circshift(trajd, 10); 
            figure(100); clf; 
            for l = 1:nqa
                subplot(nqa/2, 2, l);
                plot(trajd(:,1), trajd(:,l+1));
            end
            
            % parameters
            in_params.NSAMP = NSAMP/2.5;
            in_params.T = 1/freq(i);
            in_params.nq = nq;
            in_params.nv = nv;
            in_params.nu = nu;
            in_params.nl = nl;
            in_params.qa = qa;
            in_params.nqa = nqa;
            
            % compute inverse dynamics
            fprintf('Computing ID Solution for DC = %d, DL = %d, and F = %d Hz \r', DC(j), DL(k), freq(i)*1e3); tic;
            [xtraj, utraj, vtraj, z, F, info] = ComputeIDSol(hamrID, actuators, trajd, v, in_params);
            fprintf('Computed in %f mins \r', toc/60);
            
            %% Plotting
            ttopt = xtraj.getBreaks();
            xxopt = xtraj.eval(ttopt);
            uuopt = utraj.eval(ttopt);
            vvopt = vtraj.eval(ttopt);
            
            figure(1); clf;
            for l = 1:nqa
                subplot(nqa/2, 2, l); hold on;
                plot(ttopt, xxopt(qa(l),:))
                plot(linspace(0, in_params.T, NSAMP), trajd(:,l+1)', 'k--')
            end
            drawnow;
            
            figure(2); clf;
            NCYC = 2;
            for l = 1:nqa
                subplot(nqa/2, 2, l); hold on;
                %     plotyy(ttopt, uuopt(j,:), ttopt, vvopt(j,:))
                plot(linspace(0, NCYC, NCYC*(NSAMP/2.5)), repmat(vvopt(l,:), 1, NCYC), '*-')
                ylim([actuators.dummy_bender(1).dp.Vg, actuators.dummy_bender(1).dp.Vb])
            end
            drawnow;
%             tilefigs;
            
            %% Playback
%             ttpb = ttopt/1e3; 
%             xtrajpb = PPTrajectory(foh(ttopt, xxopt));
%             xtrajpb = xtrajpb.setOutputFrame(v.getInputFrame());
%             v.playback(xtrajpb, struct('slider', true))
%             
            %% Saving
            ffinput_jk.params.NSAMP = NSAMP;
            ffinput_jk.params.DC = DC(j);
            ffinput_jk.params.DL = DL(k);
            ffinput_jk.params.GAIT = gait;
            ffinput_jk.params.LIFTAMP = LIFTAMP;
            ffinput_jk.params.SWINGAMP = SWINGAMP;
            ffinput_jk.params.TYPE = TYPE;
            
            ffinput_jk.trajd = trajd;
            ffinput_jk.ttopt = ttopt;
            ffinput_jk.xxopt = xxopt;
            ffinput_jk.uuopt = uuopt;
            ffinput_jk.vvopt = vvopt;
            
            ffinputs{j,k} = ffinput_jk;
            
        end
        
    end  
    
    save([save_dir, 'FFInputs_Type', num2str(TYPE), '_' num2str(1e3*freq(i)), 'Hz'], 'ffinputs')
    
end

