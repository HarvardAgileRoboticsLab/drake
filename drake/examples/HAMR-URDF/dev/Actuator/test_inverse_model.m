clear; clc; close all;
addpath('./..')
load TROT_baised_225V_1Hz_100

%% Build robot
urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', 'HAMRVariational_scaled.urdf');

% options
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.z_inactive_guess_tol = 1;
options.use_bullet = false;
options.dt = 1;
options.floating = false;
options.collision = false;

options.terrain = RigidBodyFlatTerrain();
h = HamrVariationalTSRBM(urdf, options);

%% Make actuators
dp.Vb = 225;
dp.Vg = 0;
%
nact = 8;
hr_actuators = HamrActuators(nact, {'FLsact', 'FLlact', 'RLsact', 'RLlact', ...
    'FRsact', 'FRlact', 'RRsact', 'RRlact'},  [1; 1; -1; -1; 1; 1; -1; -1], dp);

% make lift's double thick
tcfL = 2*hr_actuators.dummy_bender(1).tcf; 
for i = 1:numel(hr_actuators.dummy_bender)
    if contains(hr_actuators.names{i}, 'lact')
        hr_actuators.dummy_bender(i) = hr_actuators.dummy_bender(i).setCFThickness(tcfL); 
    end
end

%%

nq = h.getNumPositions();
nu = h.getNumInputs();

xx = yy(1:2*nq,:);
uu = yy(2*nq+1:2*nq+nu,:);
vv = zeros(size(uu));

qact_v = xx(h.getActuatedJoints(),:);

for i = 1:numel(tt)
    qact = xx(h.getActuatedJoints(),i);   
%     [kact, dg] = test_dEdq(qact); disp(kact)
%     [~, dg_2] = geval(@(x)test_dEdq(x), qact, struct('grad_method', 'numerical')); 
%     valuecheck(dg, dg_2, 1e-5);
    [vv(:,i), dvv] = ipzt_fun(h, hr_actuators, tt(i), xx(:,i), uu(:,i));
end

%%
for i = 1:size(Vact,1)
    subplot(size(Vact,1)/2, 2, i)
    plot(tt, vv(i,:) - mean(vv(:,i))); hold on;
    plot(tt, Vact(i,:) - mean(Vact(:,i)), '--');
end
