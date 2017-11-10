function [xf,objval,exitflag,infeasible_constraint_name] = SLSimpleFitNLP(data, IND, SFB_INPUTS, NK, NP)

% Load dataset
N = numel(IND);                                       % Number of points in dataset
XI = data.x(:,IND);
TAU = data.tau(:,IND);

%% Build Full Single Leg

% options
name = 'FL_scaled';
SLurdf = fullfile(getDrakePath, 'examples', 'HAMR-URDF', 'dev', 'SimpleHAMR', 'urdf', [name, '.urdf']);
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = false;
options.pf = [0; 7.58; -11.35];
options.foot = 'FLL4';

SL = SLRBM(SLurdf, options);
nq = SL.getNumPositions();
nv = SL.getNumVelocities();
na = SL.getNumPositions();
nu = SL.nu; 
nl = SL.nl();

% XI = [XI(1:nq, :); rand(nv, N); XI(nq+nv+(1:nu), :)]; 

%% Build Simple Single Leg

SLSimpleurdf = fullfile(getDrakePath, 'examples', 'HAMR-URDF', 'dev', 'SimpleHAMR', ...
    'urdf','SLSimple_scaled.urdf');

% Initial guess for stiffness and damping
% init = load('SpringDamper_100Hz_80V_2_2');
K0 = -rand(numel(SFB_INPUTS), numel(SFB_INPUTS), NK);
P0 = -rand(numel(SFB_INPUTS), numel(SFB_INPUTS), NP);
% K0 = init.K;
% P0 = init.P; 
z0 = [K0(:); P0(:)];
nz = numel(z0);

% options
optionsSimple.ignore_self_collisions = true;
optionsSimple.collision_meshes = false;
optionsSimple.use_bullet = false;
optionsSimple.floating = false;
optionsSimple.k = reshape(K0, [], NK);
optionsSimple.p = reshape(P0, [], NP);

SLSimple = SLSimpleRBM(SLSimpleurdf, optionsSimple);
nqS = SLSimple.getNumPositions();
nvS = SLSimple.getNumVelocities();
naS = SLSimple.getNumPositions();
nuS = SLSimple.getNumInputs();


%% Define NLP

z_inds = reshape(1:nz, nz, 1);
a_inds = reshape(nz+ (1:(N*na)), na, N);
l_inds = reshape(nz+N*na + (1:(N*nl)), nl, N);
aS_inds = reshape(nz+N*na + N*nl + (1:(N*naS)), naS, N);

num_vars = nz + (na + naS + nl)*N;
x_name = cell(num_vars, 1);


for i = 1:nz
    x_name{i} = sprintf('z%d',i);
end

for i = 1:N
    for j = 1:na
        x_name{nz+(i-1)*na+j} = sprintf('a%d[%d]',j,i);
    end
    for j = 1:nl
        x_name{nz+N*na+(i-1)*nl+j} = sprintf('l%d[%d]',j,i);
    end
    for j = 1:naS
        x_name{nz+N*na+N*nl+(i-1)*naS+j} = sprintf('aS%d[%d]',j,i);
    end
end

nlp = NonlinearProgram(num_vars, x_name);

nlp = nlp.setSolver('snopt');
nlp = nlp.setSolverOptions('snopt','MajorIterationsLimit',10000);
nlp = nlp.setSolverOptions('snopt','MinorIterationsLimit',200000);
nlp = nlp.setSolverOptions('snopt','IterationsLimit',5000000);
nlp = nlp.setSolverOptions('snopt','SuperbasicsLimit',1000);

nlp = nlp.setSolverOptions('snopt','MajorOptimalityTolerance',1e-4);
nlp = nlp.setSolverOptions('snopt','MinorOptimalityTolerance',1e-4);
nlp = nlp.setSolverOptions('snopt','MajorFeasibilityTolerance',1e-4);
nlp = nlp.setSolverOptions('snopt','MinorFeasibilityTolerance',1e-4);
nlp = nlp.setSolverOptions('snopt','constraint_err_tol',1e-6);


%% Add Constraints and objective

copt.grad_level = 1;

% dynamics of full SL
for i = 1:N
    nvars1 = nl + na;
    dyn_constraint = FunctionHandleConstraint(zeros(na,1), zeros(na,1), nvars1, ...
        @(a, l)dyn_const_fun(a, l, XI(:,i)), copt);
    dyn_constraint = dyn_constraint.setName(sprintf('dynamic[%d]', i));
    nlp = nlp.addConstraint(dyn_constraint, {a_inds(:,i); l_inds(:,i)});
    
    % Kinematics of full SL
    acc_constraint = FunctionHandleConstraint(zeros(nl,1), zeros(nl,1), na, ...
        @(a)acc_const_fun(a, XI(1:(nq+nv),i)), copt);
    acc_constraint = acc_constraint.setName(sprintf('acceleration[%d]', i));
    nlp = nlp.addConstraint(acc_constraint, a_inds(:,i));
    
    % dynamics of simple SL
    nvars2 = naS + nz;
    dynS_constraint = FunctionHandleConstraint(zeros(naS,1), zeros(naS,1), nvars2, ...
        @(aS, z)dynS_const_fun(aS, z, [XI([SFB_INPUTS; nq+SFB_INPUTS], i); TAU(:,i)]), copt);
    dynS_constraint = dynS_constraint.setName(sprintf('dynamicS[%d]', i));
    nlp = nlp.addConstraint(dynS_constraint, {aS_inds(:,i); z_inds});
    
    % Objective
    nlp = nlp.addCost(FunctionHandleObjective(2*naS, @objective_fun), ...
        {a_inds(SFB_INPUTS, i); aS_inds(:,i)});
    
end

% NSD constraint on K and P
lb = repmat([0; -Inf], (NK+NP), 1); 
ub = repmat([Inf; 0], (NK+NP), 1); 
nsd_constraint = FunctionHandleConstraint(lb, ub, nz, @nsd_const_fun);
nsd_constraint = nsd_constraint.setName('nsd_constraint');
nlp = nlp.addConstraint(nsd_constraint, z_inds);



%% Initialize and Solve

tic;
x0 = [z0; zeros(num_vars-nz, 1)];
[xf,objval,exitflag,infeasible_constraint_name] = nlp.solve(x0);
toc

%% Objective
    function [f, df] = objective_fun(adyn, aSdyn)
        
        xin = [adyn; aSdyn];
        [f,df] = objective(xin);
        fprintf('Objective: %f \r', max(f));
        
%         step = 1e-6;
%         df_fd = zeros(size(df));
%         dxin = step*eye(length(xin));
%         for k = 1:length(xin)
%             df_fd(:,k) = (objective(xin+dxin(:,k)) - objective(xin-dxin(:,k)))/(2*step);
%         end
%         
%         disp('Objective Derivative Error:');
%         disp(max(abs(df_fd(:)-df(:))));
        
    end


    function [f, df] = objective(xin)
        
        a = xin(1:naS);
        aS = xin(naS+(1:naS));
        
        Q = eye(naS);
        f = (1/2)*(a - aS)'*Q*(a - aS);
        df = [(a - aS)'*Q, (aS - a)'*Q];
    end


    function [f, df] = dyn_const_fun(ak, lk, xk)
        xin = [ak; lk; xk];
        
        [f,df] = dyn_const(xin);
        fprintf('Full Dynamics: %f \r', max(abs(f)));
        
        %         step = 1e-6;
        %         df_fd = zeros(size(df));
        %         dxin = step*eye(length(xin));
        %         for k = 1:length(xin)
        %             df_fd(:,k) = (dyn_const(xin+dxin(:,k)) - dyn_const(xin-dxin(:,k)))/(2*step);
        %         end
        %
        %         disp('Dyn Derivative Error:');
        %         disp(max(max(abs(df_fd(:, 1:(na+nl))-df(:,1:(na+nl))))));
        
    end

    function [f, df] = dyn_const(xin)
        
        nq = SL.getNumPositions();
        nv = SL.getNumVelocities();
        nu = SL.nu;
        
        ak = xin(1:na);
        lk = xin(na+(1:nl));
        xk = xin(na+nl+(1:(nq+nv)));
        uk = xin(na+nl+nq+nv+(1:nu));
        
        
        [xkdot, dxkdot] = SL.dynamics(0, xk, [uk; lk]);
        f = xkdot(nv+(1:na)) - ak;
        df = [-eye(na), dxkdot(nv+(1:na), 1+nq+nv+nu+(1:nl))];
        
    end

    function [f, df] = dynS_const_fun(aSk, z, xSk)
        xin = [aSk; z; xSk];
        
        [f,df] = dynS_const(xin);
        fprintf('Simple Dynamics: %f \r', max(abs(f)));
        
        %         step = 1e-6;
        %         df_fd = zeros(size(df));
        %         dxin = step*eye(length(xin));
        %         for k = 1:length(xin)
        %             df_fd(:,k) = (dynS_const(xin+dxin(:,k)) - dynS_const(xin-dxin(:,k)))/(2*step);
        %         end
        %
        %         disp('Dyn Simple Derivative Error:');
        %         disp(max(max(abs(df_fd(:, 1:(naS+nz))-df(:,1:(naS+nz))))));
        %
    end


    function [f, df] = dynS_const(xin)
        
        ak = xin(1:naS);
        z = xin(naS+(1:nz));
        xk = xin(naS+nz+(1:(nqS+nvS))) + [SLSimple.q0; 0*SLSimple.q0];
        tauk = xin(naS+nz+nqS+nvS+(1:nuS));
        
        [xkdot, dxkdot_dK, dxkdot_dP] = SLSimple.dynamics_sd_fit([0; xk; tauk; z]);
        f = xkdot(nvS+(1:naS)) - ak;
        df= [-eye(naS), dxkdot_dK(nvS+(1:naS), :), dxkdot_dP(nvS+(1:naS), :)];
        
    end

    function [f, df] = acc_const_fun(ak, xk)
        
        xin = [ak; xk];
        [f,df] = acc_const(xin);
        fprintf('Accleration: %f \r', max(abs(f)));
        
        %         step = 1e-6;
        %         df_fd = zeros(size(df));
        %         dxin = step*eye(length(xin));
        %         for k = 1:length(xin)
        %             df_fd(:,k) = (acc_const(xin+dxin(:,k)) - acc_const(xin-dxin(:,k)))/(2*step);
        %         end
        %
        %         disp('Acc Const Derivative Error:');
        %         disp(max(max(abs(df_fd(:,1:na)-df(:,1:na)))));
        %
    end

    function [f, df] = acc_const(xin)
        
        ak = xin(1:na);
        qk = xin(na+(1:nq));
        vk = xin(na+nq+(1:nv));
        
        [~, dK, d2K] = SL.positionConstraints(qk);
        dK = dK(SL.valid_loops, :);
        d2K = reshape(d2K(SL.valid_loops,:)', nq, nl*nq)'; %
        f = (kron(vk', eye(nl))*d2K)*vk + dK*ak;
        df = dK;        
    end

    function [f, df] = nsd_const_fun(z)
        
        xin = z;
        [f,df] = nsd_const(xin);
%         fprintf('NSD Constraint: %f \r', max(abs(f)));
        
%         step = 1e-6;
%         df_fd = zeros(size(df));
%         dxin = step*eye(length(xin));
%         for k = 1:length(xin)
%             df_fd(:,k) = (nsd_const(xin+dxin(:,k)) - nsd_const(xin-dxin(:,k)))/(2*step);
%         end
%         
%         disp('NSD Const Derivative Error:');
%         disp(max(abs(df_fd(:)-df(:))));
        
    end

    function [f, df] = nsd_const(xin)
        k = xin(1:(NK*nqS^2)); 
        p = xin(NK*nqS^2+(1:(NP*nvS^2))); 
        
        % Constraints on K matrices
        fk = zeros(2*NK,1);
        dfk = zeros(2*NK, numel(k));        
        for ii = 1:NK
            % determinant is positive
            fk(2*(ii-1) + 1) = k(nqS^2*(ii-1) + 1)*k(nqS^2*(ii-1) + 4) - k(nqS^2*(ii-1) + 2)*k(nqS^2*(ii-1) + 3); % >=0
            dfk(2*(ii-1) + 1, nqS^2*(ii-1)+(1:nqS^2)) = [k(nqS^2*(ii-1) + 4), -k(nqS^2*(ii-1) + 3), -k(nqS^2*(ii-1) + 2), k(nqS^2*(ii-1) + 1)]; % <=0

            % trace is negative
            fk(2*(ii-1) + 2) = k(nqS^2*(ii-1) + 1) + k(nqS^2*(ii-1) + 4); % <=0
            dfk(2*(ii-1) + 2, [nqS^2*(ii-1)+1, nqS^2*(ii-1)+4]) = [1, 1]; 
        end
            
        % Constraints on P matrices
        fp = zeros(2*NP,1);
        dfp = zeros(2*NP, numel(p));        
        for ii = 1:NP
            
            % determinant is positive
            fp(2*(ii-1) + 1) = p(nvS^2*(ii-1) + 1)*p(nvS^2*(ii-1) + 4) - p(nvS^2*(ii-1) + 2)*p(nvS^2*(ii-1) + 3); % >=0
            dfp(2*(ii-1) + 1,  nvS^2*(ii-1)+(1:nvS^2)) = [p(nvS^2*(ii-1) + 4), -p(nvS^2*(ii-1) + 3), ...
                -p(nvS^2*(ii-1) + 2), p(nvS^2*(ii-1) + 1)]; % <=0

            % trace is negative
            fp(2*(ii-1) + 2) = p(nvS^2*(ii-1) + 1) + p(nvS^2*(ii-1) + 4); % <=0
            dfp(2*(ii-1) + 2, [nvS^2*(ii-1)+1, nvS^2*(ii-1)+4]) = [1, 1]; 
        end
        
        f = [fk; fp]; 
        df = blkdiag(dfk, dfp);         
        
    end

end