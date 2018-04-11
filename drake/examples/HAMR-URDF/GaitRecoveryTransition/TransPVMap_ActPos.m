function [] = TransPVMap_ActPos()

npts = 21;
lift_lim = [-0.2, 0.2];
swing_lim = [-0.2, 0.2];

% f: (FL, RL, FR, RR);
% df: (FL, RL, FR, RR);
[ss, ll, f, df] = gen_trans_data(npts, swing_lim, lift_lim);

save('ActuatorToLegPositionVelocity', 'ss', 'll', 'f', 'df')

end


function [ss, ll, f, df] = gen_trans_data(npts, slim, llim)
%% Build Full Model

urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', 'HAMR_scaledV2.urdf');

options.terrain = []; %
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = false;
options.collision = false;
options.dt = 1; %0.1;

% Build robot
hamr = HamrTSRBM(urdf, options);
nq = hamr.getNumPositions();
nv = hamr.getNumVelocities();
nc = 4; %hamr.getNumContactPairs();

all_dof = (1:nq)';
indep_dof = hamr.getActuatedJoints();
dep_dof = all_dof(sum(bsxfun(@eq, all_dof, indep_dof'),2) ~= 1);
dim = 3;


%% Leg Positions

fp_body = [0.06, 8.565, -6.322;
    -0.06, 8.565, -6.322;
    0.06, -8.565, -6.322;
    -0.06, -8.565, -6.322];

legs = {'FLL4', 'RLL4', 'FRL4', 'RRL4'};

%% Build Test Cases

% test actuator inputs
swing = linspace(slim(1), slim(2), npts);
lift = linspace(llim(1), llim(2), npts);
[ss, ll] = meshgrid(swing, lift);

% position mapping: (FL, RL, FR, RR);
f = cell(npts,npts, nc);

% velocity mapping: (FL, RL, FR, RR);
df = cell(npts,npts,nc);


%% Compute leg position and transmission ratio

yleg0 = 0*fp_body;
kinsol0 = hamr.doKinematics(zeros(nq,1), zeros(nv,1));

for k = 1:nc
    yleg0(k,:) = hamr.forwardKin(kinsol0, ...
        hamr.findLinkId(legs{k}), fp_body(k,:)')';
end

for is = 1:npts %:npts
    for il = 1:npts
        
        % initial guess
        q0 = zeros(nq, 1);
        
        %solve for depedent qs
        q_is_il = HAMR_FK(hamr, ss(is, il), ll(is, il), q0, fp_body);
        disp(q_is_il(indep_dof))
        
        % constraint gradients needed
        [~, dK] = hamr.positionConstraints(q_is_il);
        dg_dqi = dK(hamr.valid_loops, indep_dof);
        dg_dqd = dK(hamr.valid_loops, dep_dof);
        
        kinsol = hamr.doKinematics(q_is_il, 0*q_is_il);
        
        % find leg position and map from (vs, vl) to vleg (x, y, z);
        % given yleg = f(qi, qd) and g(qi, qd) = 0
        for k = 1:nc
            
            [ylegp, dylegp] = hamr.forwardKin(kinsol, ...
                hamr.findLinkId(legs{k}), fp_body(k,:)');
            
            df_dqi = dylegp(:, indep_dof);
            df_dqd = dylegp(:, dep_dof);
            Hp = df_dqi - df_dqd*inv(dg_dqd)*dg_dqi;
            
            f{is, il, k} = (ylegp' - yleg0(k,:))';
            df{is, il, k} = Hp(:, (2*k-1):2*k);
        end
    end
end

% subtract out 0 leg position
% for i = 1:nc
%     for j = 1:dim
%         yleg(:,:,j,i) = yleg(:,:,j,i) - yleg(i0, i0, j, i);
%     end
% end


%% Plotting
title_str = {'FL_X', 'FL_Y', 'FL_Z' ...
    'RL_X', 'RL_Y', 'RL_Z', ...
    'FR_X', 'FR_Y', 'FR_Z', ...
    'RR_X', 'RR_Y', 'RR_Z',};

figure(1); clf;
for i = 1:nc
    for j = 1:dim
        subplot(nc,dim, (i-1)*dim + j); hold on;
        title(title_str{(i-1)*dim + j}, 'interpreter', 'none');
        contour(ss, ll, cellfun(@(x) x(j), f(:,:,i)), 10, 'showtext', 'on');
        xlabel('Q Swing (mm)');
        ylabel('Q Lift (mm)')
    end
end


index = reshape(1:2*dim, 2, dim)';
for i = 1:nc
    figure(i+1); clf;
    for j = 1:2*dim
        subplot(dim,2,index(j)); hold on;
        [r, c] = ind2sub([dim, 2],j);
        title([legs{i}, '_', num2str(r), num2str(c)], 'interpreter', 'none');
        contour(ss, ll, cellfun(@(x) x(j), df(:,:,i)), 10, 'showtext', 'on');
        xlabel('Q Swing (mm)')
        ylabel('Q Lift (mm)')
    end
end


end

function qf = HAMR_FK(hamr, q_sact, q_lact, q0, pfLeg)

nq = hamr.getNumPositions();
nl = numel(hamr.valid_loops);
act_dof = hamr.getActuatedJoints();
nact = numel(act_dof);

q_inds = 1:nq;
num_vars = hamr.getNumPositions();
x_name = cell(num_vars, 1);

for i = 1:num_vars
    x_name{i} = sprintf('a[%d]',i);
end

nlp = NonlinearProgram(num_vars, x_name);

nlp = nlp.setSolver('snopt');
nlp = nlp.setSolverOptions('snopt','MajorIterationsLimit',10000);
nlp = nlp.setSolverOptions('snopt','MinorIterationsLimit',200000);
nlp = nlp.setSolverOptions('snopt','IterationsLimit',5000000);
nlp = nlp.setSolverOptions('snopt','SuperbasicsLimit',1000);

tol = 1e-6;
nlp = nlp.setSolverOptions('snopt','MajorOptimalityTolerance',tol);
nlp = nlp.setSolverOptions('snopt','MinorOptimalityTolerance',tol);
nlp = nlp.setSolverOptions('snopt','MajorFeasibilityTolerance',tol);
nlp = nlp.setSolverOptions('snopt','MinorFeasibilityTolerance',tol);
nlp = nlp.setSolverOptions('snopt','constraint_err_tol',tol);

% Loop Constraints
loop_const = FunctionHandleConstraint(zeros(nl, 1), zeros(nl, 1), num_vars, ...
    @loop_constraint);
loop_const = loop_const.setName('loop_constraint');
nlp = nlp.addConstraint(loop_const, q_inds);

% Set Leg Position
des_act_pos = (ones(nact/2, 1)*[q_sact, q_lact])';
nlp = nlp.addConstraint(ConstantConstraint(des_act_pos), q_inds(act_dof));

% Bounding box Constraints
[q_min, q_max] = getJointLimits(hamr);      % joints limits
nlp = nlp.addConstraint(BoundingBoxConstraint(1.1*q_min, 1.1*q_max), q_inds);



% solve
qf = nlp.solve(q0);

    function [f, df] = loop_constraint(q)
        manip = hamr.getManipulator();
        [k, dK] = manip.positionConstraints(q);
        f = k(hamr.valid_loops);
        df = dK(hamr.valid_loops,:);
    end
end
