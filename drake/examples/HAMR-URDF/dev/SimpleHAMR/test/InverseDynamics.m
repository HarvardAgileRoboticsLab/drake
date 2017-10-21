function [xf,objval,exitflag,infeasible_constraint_name] = inverse_dynamics(x0, qf, vf)
%% Build Single Leg

sl_urdf = fullfile(getDrakePath, 'examples', 'HAMR-URDF', 'dev', 'SimpleHAMR', 'SLSimple_scaled.urdf');

% options
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = false; 

SL = RigidBodyManipulator(sl_urdf, options); 
SL = SL.setGravity([0; 0; -9.81e-3]); 
SL = compile(SL); 


nq = SL.getNumPositions();
nv = SL.getNumVelocities();
na = nv; 
nu = SL.getNumInputs();
nl = SL.getNumStateConstraints();

%% Define NLP

num_vars = nq + nv + na + nu + nl;
x_name = cell(num_vars, 1);
for i = 1:num_vars
    if i <= nq
        x_name{i} = sprintf('x[%d]',i);
    elseif i > nq && i <= (nq + nv)
        x_name{i} = sprintf('v[%d]',i-nq);
    elseif i > nq + nv && i <= (nq + nv + na)
        x_name{i} = sprintf('a[%d]',i-nq-nv);    
    elseif i > (nq+nv+na) && i <= (nq + nv + na + nu)
        x_name{i} = sprintf('u[%d]',i-nq-nv-na);
    elseif i > (nq+nv+nu +na) && i <= (nq + nv + na + nu +nl)
        x_name{i} = sprintf('l[%d]',i-nq-nv-na-nu);
    end
end


nlp = NonlinearProgram(num_vars, x_name);

% objective
nlp = nlp.addCost(FunctionHandleObjective(nq+nv,@(x)objective_fun(SL,x)), (1:(nq+nv))');

% input constraint
ulim = 10;
nlp = nlp.addConstraint(BoundingBoxConstraint(-ulim*ones(nu,1), ulim*ones(nu,1)), (nq+nv+na+(1:nu))');

% state constraint
[qmin, qmax] = SL.getJointLimits();
nlp = nlp.addConstraint(BoundingBoxConstraint(qmin, qmax), (1:nq)');

% dynamics constraint
dyn_constraint = FunctionHandleConstraint(zeros(nq,1), zeros(nq,1), num_vars, @(x)dyn_const_fun(SL, x)); 
dyn_constraint = dyn_constraint.setName('dynamics'); 
nlp = nlp.addConstraint(dyn_constraint, (1:num_vars)'); 


%% Initialize and Solve

tic;
[xf,objval,exitflag,infeasible_constraint_name] = nlp.solve(x0);
toc


%% Const and Objectives

    function [f, df] = objective_fun(obj, xin)
        
        [f,df] = objective(obj,xin);        
        
        df_fd = zeros(size(df));
        step = sqrt(eps(max(xin)));
        dxin = step*eye(length(xin));
        for k = 1:length(xin)
            df_fd(:,k) = (objective(obj, xin+dxin(:,k)) - objective(obj, xin-dxin(:,k)))/(2*step);
        end        
        
        disp('Objective Derivative Error:');
        disp(max(abs(df_fd(:)-df(:))));
        
    end

    function [f, df] = objective(obj, x)
        
        qk = x(1:nq); 
        vk = x(nq+(1:nv)); 
        
        pf = [0 0 -14.97]'; % position of foot in local frame
        kinsolk = obj.doKinematics(qk, vk, struct('compute_gradients', true));        
        [xfoot, Jk, dJk] = obj.forwardKin(kinsolk, obj.findLinkId('FL2'), pf);   
        vfoot = Jk*qk;
                
        f1 = (1/2)*(xfoot - qf)'*(xfoot - qf);      % position 
%         f2 = 0; 
        f2 = (1/2)*(vfoot - vf)'*(vfoot - vf);      % velocity
        f2 = 0; 
        f = f1 + f2; 
        
               
        df1 = [(xfoot - qf)'*Jk, zeros(1,nv)];      % df/dq, df/dv
        df2 = [kron(vk', (vfoot - vf)')*reshape(dJk, 3*nv, nv), (vfoot - vf)'*Jk]; 
        df2 = 0; 
        df = df1 + df2; 
        
    end

    function [f, df] = dyn_const_fun(obj, xin)
        
        [f,df] = dyn_const(obj,xin);        
        
%         df_fd = zeros(size(df));
%         step = sqrt(eps(max(xin)));
%         dxin = step*eye(length(xin));
%         for k = 1:length(xin)
%             df_fd(:,k) = (dyn_const(obj, xin+dxin(:,k)) - dyn_const(obj, xin-dxin(:,k)))/(2*step);
%         end
%         
%         
%         disp('Dynamics Derivative Error:');
%         disp(max(abs(df_fd(:)-df(:))));
    end

    function [f, df] = dyn_const(obj, x)
        
        qk = x(1:nq); 
        vk = x(nq+(1:nv)); 
        ak = x(nq+nv+(1:na)); 
        uk = x(nq+nv+na+(1:nu)); 
%         
%         kinsolk = obj.doKinematics(qk, vk, struct('compute_gradients', true)); 
%         [~, Jk, dJk] = obj.forwardKin(kinsolk, obj.findLinkId('FL2'), pf);

        [Hk,Ck,Bk,dHk,dCk,dBk] = obj.manipulatorDynamics(qk, vk);
        
        f = Hk*ak + Ck - Bk*uk; 
        df = [kron(ak', eye(nq))*dHk(:,1:nq) + dCk(:,1:nq) - kron(uk', eye(nq))*dBk(:,1:nq), ...  % df/dq 
            dCk(:,nq+(1:nv)) -  kron(uk', eye(nv))*dBk(:,nq+(1:nv)), ...                       % df/dv
            Hk, ...                                                        % df/da
            -Bk];                                                      % df/du                        
    end


end

