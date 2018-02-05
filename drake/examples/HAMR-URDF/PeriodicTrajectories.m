clear; clc; close all;
%% Load Files

load_dir = '~/Dropbox/CurrentWork/FrictionTrajOpt/MatFiles/SimWarmStart/';
save_dir = '~/Dropbox/CurrentWork/FrictionTrajOpt/MatFiles/forXPC/';


% files = {'TROT_0.15N_30Hz_TYM_SP_VariationalSmooth'}; 

% files = {'TROT_0.15N_30Hz_TYM_CS_VariationalSmooth', 'TROT_0.25N_30Hz_TYM_SP_VariationalSmooth'};
files = {'TROT_0.1N_30Hz_TYM_TEF_VariationalSmooth_converted'};

dim = 3;

%% Build Robot
urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', 'HAMR_scaledV2_TYM.urdf');

options.terrain = RigidBodyFlatTerrain();
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = true; %false;
options.collision = true; %false;
options.dt = 1; %0.1;

% Build robot + visualizer
hamr = HamrTSRBM(urdf, options);
nq = hamr.getNumPositions();
nv = hamr.getNumVelocities();
nu = hamr.getNumInputs();
nc = hamr.getNumContactPairs();

%% Foot and marker positions

pfm = [0.06, 8.565, -6.322;
    -0.06, 8.565, -6.322;
    0.06, -8.565, -6.322;
    -0.06, -8.565, -6.322];

pfl = [0, 7.540, -11.350;
    0, 7.540, -11.350;
    0, -7.540, -11.350;
    0, -7.540, -11.350];

legs = {'FLL4', 'RLL4', 'FRL4', 'RRL4'};
fkopt.base_or_frame_id = hamr.findLinkId('Chassis');

%%

title_str = {'Com-X(mm)', 'Com-Y(mm)', 'Com-Z(mm)', 'Roll(deg)', 'Pitch(deg)', 'Yaw(deg)'};

figure(1); clf;
for k = 1:numel(files)
    
    % Unpack trajectories
    trajk = load([load_dir, files{k}, '_PlusAct.mat']);
    tt = 1e-3*trajk.xtraj.getBreaks();
    hh = mean(diff(tt))*1e3;
    xx = trajk.xtraj.eval(tt*1e3);
    vv = trajk.xtraj.eval(tt*1e3 + hh/2);
    uu = trajk.vtraj.eval(tt*1e3 + hh/2);
    ttv = tt + hh*1e-3/2;
    
    % Plot floating base configuration
    figure(1); hold on;
    for j = 1:2*dim
        subplot(2,dim,j); hold on; ylabel(title_str(j), 'FontSize', 18)
        xlabel('Time(s)', 'FontSize', 18)
        if j > 3
            plot(tt, rad2deg(xx(j,:)), 'LineWidth', 1.5);
        else
            plot(tt, xx(j,:), 'LineWidth', 1.5);
        end
    end
    
    % Compute leg and marker positions
    pf0 = 0*pfm;
    pfW = zeros([numel(tt), size(pfm')]);
    pfW2 = zeros([numel(tt), size(pfl')]);
    pfB = zeros([numel(tt), size(pfm')]);
    pfB2 = zeros([numel(tt), size(pfl')]);
    vfB = zeros([numel(tt), size(pfm')]);
    
    
    kinsol0 = hamr.doKinematics(zeros(nq,1), zeros(nv,1));
    for i = 1:size(pfm,1)
        pf0(i,:) = hamr.forwardKin(kinsol0, hamr.findLinkId(legs{i}), pfm(i,:)')';
    end
    
    for i = 1:size(pfm,1)
        for j = 1:numel(tt)
            q = xx(1:nq, j);
            qd = xx(nq+(1:nv), j);
            kinsol = hamr.doKinematics(q, qd, struct('compute_gradients', true));
            
            pfW(j,:,i) = hamr.forwardKin(kinsol, hamr.findLinkId(legs{i}), ...
                pfm(i,:)');            
            [pfB(j,:,i), dpfB_dx] = hamr.forwardKin(kinsol, hamr.findLinkId(legs{i}), ...
                pfm(i,:)',fkopt);                        
            pfB(j,:,i) = pfB(j,:,i) - pf0(i,:);
            vfB(j,:,i) = dpfB_dx*qd;
            
            pfW2(j,:,i) = hamr.forwardKin(kinsol, hamr.findLinkId(legs{i}), ...
                pfl(i,:)');
            [pfB2(j,:,i), ~] = hamr.forwardKin(kinsol, hamr.findLinkId(legs{i}), ...
                pfl(i,:)',fkopt);

        end
    end
    
    figure(2); hold on;
    for j = 1:nc
        
        subplot(nc,dim,dim*(j-1)+1); hold on;
        ylabel('Foot X'); title(legs{j})
        plot(repmat(pfB(:,1,j), 2, 1),'LineWidth', 1.5)
        
        subplot(nc,dim,dim*(j-1)+2); hold on;
        ylabel('Foot Y'); title(legs{j})
        plot(repmat(pfB(:,2,j), 2, 1),'LineWidth', 1.5)
        
        subplot(nc,dim,dim*(j-1)+3); hold on;
        ylabel('Foot Z'); title(legs{j})
        plot(repmat(pfB(:,3,j), 2, 1),'LineWidth', 1.5)
        
    end
    
    figure(3); hold on;
    for j = 1:nc
        
        subplot(nc,dim,dim*(j-1)+1); hold on;
        ylabel('Foot X'); title(legs{j})
        plot(repmat(pfW(:,1,j), 2, 1),'LineWidth', 1.5)
        
        subplot(nc,dim,dim*(j-1)+2); hold on;
        ylabel('Foot Y'); title(legs{j})
        plot(repmat(pfW(:,2,j), 2, 1),'LineWidth', 1.5)
        
        subplot(nc,dim,dim*(j-1)+3); hold on;
        ylabel('Foot Z'); title(legs{j})
        plot(repmat(pfW(:,3,j), 2, 1),'LineWidth', 1.5)
        
    end
    
    figure(4); hold on;
    for j = 1:nc
        
        subplot(nc,dim,dim*(j-1)+1); hold on;
        ylabel('Foot Vel X'); title(legs{j})
        plot(repmat(vfB(:,1,j),2,1),'LineWidth', 1.5)
        
        subplot(nc,dim,dim*(j-1)+2); hold on;
        ylabel('Foot Vel Y'); title(legs{j})
        plot(repmat(vfB(:,2,j),2,1),'LineWidth', 1.5)
        
        subplot(nc,dim,dim*(j-1)+3); hold on;
        ylabel('Foot Vel Z'); title(legs{j})
        plot(repmat(vfB(:,3,j),2,1),'LineWidth', 1.5)
        
    end
    
    save([save_dir, files{k}, '_TYM_forXPC.mat'], 'tt', 'ttv', 'xx', 'uu', ...
        'pfW',  'pfB', 'vfB', 'pfW2', 'pfB2');
%     
end


