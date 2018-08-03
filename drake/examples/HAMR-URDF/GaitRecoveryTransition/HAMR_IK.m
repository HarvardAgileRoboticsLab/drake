function qf = HAMR_IK(hamr, legx, legz, q0)

nq = hamr.getNumPositions();
nl = numel(hamr.VALID_LOOPS);
act_dof = hamr.getActuatedJoints();
nact = numel(act_dof);

q_inds = 1:nq;
num_vars = hamr.getNumPositions();
x_name = cell(num_vars, 1);

for i = 1:num_vars
    x_name{i} = sprintf('q[%d]',i);
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
fpopt.base_frame = 'World';
fpopt.loc = 'foot';
pf0 = hamr.getFootPosition(0*q0, 0*q0, fpopt);

% pf0 = hamr.ge
% pf0 = 0*pfLeg;
% legs = hamr.legs;
% kinsol0 = hamr.doKinematics(zeros(nq,1), zeros(nq,1));
% for i = 1:size(pfLeg,1)
%     pf0(i, :) = hamr.forwardKin(kinsol0, hamr.findLinkId(legs{i}),...
%         pfLeg(i,:)');
% end

% solve
qf = nlp.solve(q0);


    function [f, df] = loop_constraint(q)
        [k, dK] = hamr.positionConstraints(q);
        f = k(hamr.VALID_LOOPS);
        df = dK(hamr.VALID_LOOPS,:);
    end

    function [f, df] = leg_pos_constraint(q)
        
        f = zeros(nact, 1);
        df = zeros(nact, nq);       
        
        
        [pfk, dpfK] = hamr.getFootPosition(q, 0*q, fpopt);

        for j = 1:length(hamr.LEG_NAME)
            f((2*j-1):(2*j)) = [pfk(1,j) - pf0(1,j) - legx; ...
                pfk(3,j) - pf0(3,j) - legz];
            df((2*j-1):(2*j),:) = dpfK{j}([1;3], :);
        end
    end
end
