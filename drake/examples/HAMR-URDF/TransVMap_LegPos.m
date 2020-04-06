function [] = TransVMap_LegPos()

npts = 21;
xlim = [-3, 3];
zlim = [-1.5, 3];

% yleg: 3rd dim (x, y, z) and 4th dim (FL, RL, FR, RR);
% H: 3rd dim (FL, RL, FR, RR);
[xx, zz, qact, H] = gen_trans_data(npts, xlim, zlim);

% save('ActuatorToLegVelocity', 'xx', 'zz', 'H')

end


function [xx, zz, qact, H] = gen_trans_data(npts, xlim, zlim)
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
hamr = HamrTSRBM(HamrRBM(urdf, options), options.dt, options);
nq = hamr.getNumPositions();
nv = hamr.getNumVelocities();
nc = 4; %hamr.getNumContactPairs();

all_dof = (1:nq)';
indep_dof = hamr.getActuatedJoints();
dep_dof = all_dof(sum(bsxfun(@eq, all_dof, indep_dof'),2) ~= 1);
nact = numel(indep_dof);
dim = 3;


%% Leg Positions

% fp_body = [0.06, 8.565, -6.322;
%     -0.06, 8.565, -6.322;
%     0.06, -8.565, -6.322;
%     -0.06, -8.565, -6.322];
fp_body = [0, 7.58, -11.350;
    0, 7.58, -11.350;
    0, -7.58, -11.350;
    0, -7.58, -11.350];

legs = {'FLL4', 'RLL4', 'FRL4', 'RRL4'};

%% Build Test Cases

% test actuator inputs
x = linspace(xlim(1), xlim(2), npts);
z = linspace(zlim(1), zlim(2), npts);
[xx, zz] = meshgrid(x, z);

% actuator positons: 3rd dim (FL_swing, FL_lift, RL, FR, RR)
qact = zeros(npts,npts, nact);

% velocity mapping: 3rd dim (x, y, z) and 4th dim (FL, RL, FR, RR);
H = cell(npts,npts,nc);


%% Compute leg position and transmission ratio

% initial guess
q0 = zeros(nq, 1);
i0 = (npts-1)/2 + 1;        % index of 0,0

yleg0 = 0*fp_body;
kinsol0 = hamr.doKinematics(zeros(nq,1), zeros(nv,1));

for k = 1:nc
    yleg0(k,:) = hamr.forwardKin(kinsol0, ...
        hamr.findLinkId(legs{k}), fp_body(k,:)')';
end

for is = 1:npts %:npts
    for il = 1:npts
        
        %solve for depedent qs
        q_is_il = HAMR_IK(hamr, xx(is, il), zz(is, il), q0, fp_body);
        %         v.inspector([q_is_il; 0*q_is_il]);
        
        % constraint gradients needed
        [~, dK] = hamr.positionConstraints(q_is_il);
        dg_dqi = dK(hamr.HamrRBM.VALID_LOOPS, indep_dof);
        dg_dqd = dK(hamr.HamrRBM.VALID_LOOPS, dep_dof);
        
        
        qact(is, il, :) = q_is_il(indep_dof);
        kinsol = hamr.doKinematics(q_is_il, 0*q_is_il);
        
        % find leg position and map from (vs, vl) to vleg (x, y, z);
        % given yleg = f(qi, qd) and g(qi, qd) = 0
        for k = 1:nc
            
            [ylegp, dylegp] = hamr.forwardKin(kinsol, ...
                hamr.findLinkId(legs{k}), fp_body(k,:)');
            disp(ylegp' - yleg0(k,:))
            
            df_dqi = dylegp(:, indep_dof);
            df_dqd = dylegp(:, dep_dof);
            Hp = df_dqi - df_dqd*inv(dg_dqd)*dg_dqi;
            H{is, il, k} = Hp(:, (2*k-1):2*k);
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
title1 = {'FL_Swing', 'FL_Lift' ...
    'RL_Swing', 'RL_Lift', ...
    'FR_Swing', 'FR_Lift', ...
    'RR_Swing', 'RR_Lift',};

figure(1); clf;
for j = 1:nact
    subplot(2,nact/2,j); hold on;
    title(title1{j}, 'interpreter', 'none');
    contour(xx, zz, qact(:,:,j), 10, 'showtext', 'on');
    xlabel('Leg X (mm)');
    ylabel('Leg Z (mm)')
end


index = reshape(1:2*dim, 2, dim)';
for i = 1:nc
    figure(i+1); clf;
    for j = 1:2*dim
        subplot(dim,2,index(j)); hold on;
        [r, c] = ind2sub([dim, 2],j);
        title([legs{i}, '_', num2str(r), num2str(c)], 'interpreter', 'none');
        contour(xx, zz, cellfun(@(x) x(j), H(:,:,i)), 10, 'showtext', 'on');
        xlabel('Leg X (mm)')
        ylabel('Leg Z (mm)')
    end
end


end

function qf = HAMR_IK(hamr, legx, legz, q0, pfLeg)

nq = hamr.getNumPositions();
nl = numel(hamr.HamrRBM.VALID_LOOPS);
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
leg_pos_const = FunctionHandleConstraint(zeros(nact, 1), zeros(nact, 1), num_vars, ...
    @leg_pos_constraint);
leg_pos_const = leg_pos_const.setName('leg_pos_constraint');
nlp = nlp.addConstraint(leg_pos_const, q_inds);

% Bounding box Constraints
[q_min, q_max] = getJointLimits(hamr);      % joints limits
nlp = nlp.addConstraint(BoundingBoxConstraint(1.1*q_min, 1.1*q_max), q_inds);

% Leg initial positions
pf0 = 0*pfLeg;
legs = hamr.HamrRBM.LEG_NAME;
kinsol0 = hamr.doKinematics(zeros(nq,1), zeros(nq,1));
for i = 1:size(pfLeg,1)
    pf0(i, :) = hamr.forwardKin(kinsol0, hamr.findLinkId(legs{i}),...
        pfLeg(i,:)');
end

% solve
qf = nlp.solve(q0);


    function [f, df] = loop_constraint(q)
        manip = hamr.getManipulator();
        [k, dK] = manip.positionConstraints(q);
        f = k(hamr.HamrRBM.VALID_LOOPS);
        df = dK(hamr.HamrRBM.VALID_LOOPS,:);
    end

    function [f, df] = leg_pos_constraint(q)
        
        f = zeros(nact, 1);
        df = zeros(nact, nq);
        
        kinsol = hamr.doKinematics(q, 0*q);
        for j = 1:size(pfLeg,1)
            [pfk, dpfK] = hamr.forwardKin(kinsol, hamr.findLinkId(legs{j}),...
                pfLeg(j,:)');
            f((2*j-1):(2*j)) = [pfk(1) - pf0(j,1) - legx; ...
                pfk(3) - pf0(j,3) - legz];
            df((2*j-1):(2*j),:) = dpfK([1;3], :);
        end
        %         max(abs(f))
    end
end
