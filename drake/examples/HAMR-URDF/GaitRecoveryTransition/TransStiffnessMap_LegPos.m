function K = TransStiffnessMap_LegPos(opt)
close all;
clc;
addpath('../')


npts = 15;
xx = linspace(-2, 2, npts);
zz = linspace(-1.125, 0.375, npts);

%% Load Rigid Body Model

urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', 'HAMR_scaledV2.urdf');

% options
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.z_inactive_guess_tol = 0.1;
options.use_bullet = false;

% options to change
dt = 0.2;
options.dt = dt;
options.floating = false;
options.collision = false;

% Build robot + visualizer
hamr = HamrRBM(urdf, options);
nq = hamr.getNumPositions();
nqi = hamr.getNumInputs();


%% turn of gravity
hamr = hamr.setGravity(zeros(3,1));
hamr = compile(hamr);

%% Build actuators
dp.Vb = 225;
dp.Vg = 0;

hr_actuators = HamrActuators(nqi, {'FLsact', 'FLlact', 'RLsact', 'RLlact', ...
    'FRsact', 'FRlact', 'RRsact', 'RRlact'}, [], dp);


%% Build stiffness map
K = cell(npts,npts);

for i = 1:npts
    for j = 1:npts
        K{j,i} = gen_stiffness_data(hamr, hr_actuators, xx(i), zz(j), opt);
    end
end
save(['TransmissionStiffness_', opt{1} '_', opt{2}], 'xx', 'zz', 'K')
%% Plotting

title_str = {'FL_X', 'FL_Y', 'FL_Z' ...
    'RL_X', 'RL_Y', 'RL_Z', ...
    'FR_X', 'FR_Y', 'FR_Z', ...
    'RR_X', 'RR_Y', 'RR_Z',};

figure(100); clf; hold on;
plot(zz, cellfun(@(x) x(2,2,1), K(3,:,:)))
plot(zz, cellfun(@(x) x(2,2,1), K(:,3,:)))


for i = 1:4
    figure(i); clf;
    subplot(1,2,1); hold on; title(title_str{(i-1)*3+1})
    contourf(xx, zz, cellfun(@(x) x(1,1,i), K(:,:,:)), 10);
    %     contourf(cellfun(@(x) x(1,1,i), temp(:,:,:)), ...
    %         cellfun(@(x) x(2,1,1), temp(:,:,:)), ...
    %         cellfun(@(x) x(1,1,i), K(:,:,:)), 10);
    xlabel('Leg X (mm)')
    ylabel('Leg Z (mm)')
    colorbar;
    
    subplot(1,2,2); hold on; title(title_str{(i-1)*3+3})
    contourf(xx, zz, cellfun(@(x) x(2,2,i), K(:,:,:)), 10);
    xlabel('Leg X (mm)')
    colorbar;
end
tilefigs;

end


function K = gen_stiffness_data(hamr, hr_actuators, legx, legz, opt)

global Kp

legp = [legx; legz];
nqi = length(legp);
nact = numel(hamr.getActuatedJoints());

step = 1e-5;
dkm = zeros(nact, nqi);
dkc = zeros(nact, nqi);
dleg = step*eye(nqi);

switch opt{1}
    case 'OL'
        disp('Open Loop Stiffness')
        Kp = zeros(nact); 
    case 'CL'        
        AvgLQRGains = load('AverageLQRGains.mat');  
        
        switch opt{2}
            case '25'
                disp('Min Closed Loop Stiffness')
                Kp = AvgLQRGains.L_m1(:, 1:2:end);                
            case '50'                
                disp('Avg Closed Loop Stiffness')
                Kp = AvgLQRGains.L_ave(:, 1:2:end);                
            case '75'                
                disp('Max Closed Loop Stiffness')
                Kp = AvgLQRGains.L_p1(:, 1:2:end);
        end      
        Kp(abs(Kp) < 1e-6) = 0; 
end

for i = 1:nqi
    [Gp, Fp, ~] = spring_force_cart(hamr, hr_actuators, legp+dleg(:,i));
    [Gm, Fm, ~] = spring_force_cart(hamr, hr_actuators, legp-dleg(:,i));
    dkm(:,i) = (Gp - Gm)/(2*step);
    dkc(:,i) = (Fp - Fm)/(2*step); 
end

[~, ~, Jqi_xz] = spring_force_cart(hamr, hr_actuators, legp);
% KpF = dFdV*Kp; 

K = zeros(nqi, nqi,  nact/nqi);
for i = 1:nact/nqi
    switch opt{1}
        case 'OL'
            Ki = dkm(nqi*(i-1)+1:nqi*i, :);
        case 'CL'
            Kci = dkc(nqi*(i-1)+1:nqi*i, :);
            Kmi = dkm(nqi*(i-1)+1:nqi*i, :);
            Ki =  Kmi + Kci;
    end
    lastwarn('')
    K(:,:,i) = 1e3*((Jqi_xz{i}')\Ki);
    [warnMsg, ~] = lastwarn;
    if ~isempty(warnMsg)
        disp(legp)
        disp(Ki)
    end
end

end

function [Gqi, F, Jqi_xz] = spring_force_cart(hamr, hr_actuators, legp)

global Kp

legx = legp(1);
legz = legp(2);

nq = hamr.getNumPositions();
q0 = zeros(nq, 1);

all_dof = (1:nq)';
indep_dof = hamr.getActuatedJoints();
dep_dof = all_dof(sum(bsxfun(@eq, all_dof, indep_dof'),2) ~= 1);

nqi = length(indep_dof); 

%% find current state
qi = HAMR_IK(hamr, legx, legz, q0(1:nq));
[~,G] = hamr.manipulatorDynamics(qi, 0*qi);                          % should be G since v = 0

% constraint gradients needed
[~, K] = hamr.positionConstraints(qi);
Kqi = K(hamr.VALID_LOOPS, indep_dof);
Kqd = K(hamr.VALID_LOOPS, dep_dof);
Jqd_qi = -(Kqd\Kqi);

% Force to Voltage Mapping 
F = hr_actuators.output(0, [], [hr_actuators.dummy_bender(1).dp.Vb/2*ones(nqi,1)+Kp*qi(indep_dof); qi(indep_dof)]);
%dFdV = dFdx(:, 1+(1:nqi));
%F

% mapping from independent generalized coordinates to xyz
fpopt.base_frame = 'World';
fpopt.loc = 'foot';

% pleg0 = hamr.getFootPosition(0*qi, 0*qi, fpopt);
[~, Jleg] = hamr.getFootPosition(qi, 0*qi, fpopt);

Jqi_xz = cell(size(Jleg));
for i = 1:numel(Jleg)
    Jqi_xyz = Jleg{i}(:,indep_dof) + Jleg{i}(:,dep_dof)*Jqd_qi;
    Jqi_xz{i} = Jqi_xyz([1,3], (i-1)*length(legp)+1:i*length(legp));
    %     leg_xz(:,:,i) = pleg([1,3], i) - pleg0([1,3], i);
end
Gqi = G(indep_dof) + Jqd_qi'*G(dep_dof);

end




