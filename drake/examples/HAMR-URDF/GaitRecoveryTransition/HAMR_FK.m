function qf = HAMR_FK(hamr, q_sact, q_lact, q0)

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

% Set act pos
des_act_pos = (ones(nact/2, 1)*[q_sact, q_lact])';
nlp = nlp.addConstraint(ConstantConstraint(des_act_pos(:)'), q_inds(act_dof));

% Bounding box Constraints
[q_min, q_max] = getJointLimits(hamr);      % joints limits
nlp = nlp.addConstraint(BoundingBoxConstraint(1.1*q_min, 1.1*q_max), q_inds);


% solve
qf = nlp.solve(q0);

    function [f, df] = loop_constraint(q)
%         if is 
        if isa(hamr, 'HamrTSRBM')
            manip = hamr.getManipulator();
        else
            manip = hamr;
        end
        [k, dK] = manip.positionConstraints(q);
        f = k(hamr.VALID_LOOPS);
        df = dK(hamr.VALID_LOOPS,:);
    end
end
