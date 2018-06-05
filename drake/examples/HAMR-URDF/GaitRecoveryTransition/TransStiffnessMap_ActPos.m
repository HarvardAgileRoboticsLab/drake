function [] = TransStiffnessMap_ActPos()
close all;

npts = 11;
ss = linspace(-0.2, 0.2, npts);
ll = linspace(-0.2, 0.2, npts);

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


%% turn of gravity
hamr = hamr.setGravity(zeros(3,1));
hamr = compile(hamr);

%% Build stiffness map
K = cell(npts,npts);

for i = 1:npts
    for j = 1:npts
        K{j,i} = gen_stiffness_data(hamr, ss(i), ll(j));
    end
end


%% Plotting

title_str = {'FL_X', 'FL_Y', 'FL_Z' ...
    'RL_X', 'RL_Y', 'RL_Z', ...
    'FR_X', 'FR_Y', 'FR_Z', ...
    'RR_X', 'RR_Y', 'RR_Z',};

for i = 1:4
    figure(i); clf;
    subplot(1,3,1); hold on; title(title_str{(i-1)*3+1})
    contourf(ss, ll, cellfun(@(x) x(1,1,i), K(:,:,:)), 10);
    ylabel('K Lift (mm)')
    colorbar;
    
    
    subplot(1,3,2); hold on; title(title_str{(i-1)*3+2})
    contourf(ss, ll, cellfun(@(x) x(2,2,i), K(:,:,:)), 10);
    xlabel('K Swing (mm)');
    colorbar;
    
    
    subplot(1,3,3); hold on; title(title_str{(i-1)*3+3})
    contourf(ss, ll, cellfun(@(x) x(3,3,i), K(:,:,:)), 10);
    colorbar;
end
tilefigs;
%% Debugging

v = hamr.constructVisualizer();
v.inspector(zeros(hamr.getNumContStates(),1))

end


function K = gen_stiffness_data(hamr, q_sact, q_lact)

nq = hamr.getNumPositions();
dim = 3;
q0 = zeros(nq, 1);

all_dof = (1:nq)';
indep_dof = hamr.getActuatedJoints();
dep_dof = all_dof(sum(bsxfun(@eq, all_dof, indep_dof'),2) ~= 1);

%% leg geometry
lp_b = [0, 7.58, -11.350;
    0, 7.58, -11.350;
    0, -7.58, -11.350;
    0, -7.58, -11.350];
legs = {'FLL4', 'RLL4', 'FRL4', 'RRL4'};
nc = size(lp_b, 1);

%% find current state
qi = HAMR_FK(hamr, q_sact, q_lact, q0(1:nq), NaN);

[~,G,~,~,dG,~] = hamr.manipulatorDynamics(qi, 0*qi);                          % should be G since v = 0

% constraint gradients needed
[~, K, dK] = hamr.positionConstraints(qi);
Kqi = K(hamr.valid_loops, indep_dof);
Kqd = K(hamr.valid_loops, dep_dof);
dqi_dqd = -Kqd\Kqi;


d2V_dqi2 = dG(indep_dof, indep_dof);
d2V_dqd2 = dG(dep_dof, dep_dof);
d2V_dqdqi  = dG(dep_dof, indep_dof);


%G_qi = G(indep_dof) + (G(dep_dof)'*dqi_dqd)'

K = zeros(dim, dim, nc);
for i = 1:nc
    kinsol = hamr.doKinematics(qi, 0*qi, struct('compute_gradients', true));
    [f, df, d2f] = hamr.forwardKin(kinsol, ...
        hamr.findLinkId(legs{i}), lp_b(i,:)');
    
    df_dqi = df(:, indep_dof);
    df_dqd = df(:, dep_dof);
    
    K = d2V_dqi2 + d2V_dqdqi*dqi_dqd
    
%     df_qi = df_dqi - df_dqd*inv(Kqd)*Kqi;
%     G_xyz = df_qi*G_qi
%     K(:,:,i) = 1e3*inv(df*inv(dG_dq)*df');
    %1e3*((df*(dG_dq\df'))\eye(3));                                             % stiffness in N/m
    %     disp(K(:,:,i));
    %     disp()
end


end




