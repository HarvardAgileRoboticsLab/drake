function [xf, objval, SL] = SLSimpleIK(data, IND, SFB_INDS)

% Load dataset
N = numel(IND);                                       % Number of points in dataset
XI = data.x(:,IND);

%% Build Full Single Leg

% options
name = 'FL_scaled';
SLurdf = fullfile(getDrakePath, 'examples', 'HAMR-URDF',  'urdf', [name, '.urdf']);
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = false;
options.pf = [0; 7.58; -11.35];
options.foot = 'FLL4';

SL = SLRBM(SLurdf, options);
nq = SL.getNumPositions();
nv = SL.getNumVelocities();

%% Build foot trajectory
xyz_foot = zeros(3, N);
vxyz_foot = zeros(3,N); 
for i = 1:N
    q = XI(1:nq, i);
    qd = XI(nq+(1:nv), i);
    kinsol = SL.doKinematics(q, qd, struct('compute_gradients', true));
    [xyz_foot(:,i), Jfoot] = SL.forwardKin(kinsol, SL.findLinkId(options.foot), options.pf);
    vxyz_foot(:,i) = Jfoot*qd; 
end

%% Build Simple Single Leg

SLSimpleurdf = fullfile(getDrakePath, 'examples', 'HAMR-URDF','urdf','SLSimple_scaled.urdf');

% Initial guess for stiffness and damping
K0 = repmat(-rand(1)*eye(2), 1, 1, 1);
P0 = repmat(-rand(1)*eye(2), 1, 1, 1);

% options
optionsSimple.ignore_self_collisions = true;
optionsSimple.collision_meshes = false;
optionsSimple.use_bullet = false;
optionsSimple.floating = false;
optionsSimple.k = reshape(K0, [], 1);
optionsSimple.p = reshape(P0, [], 1);

SLSimple = SLSimpleRBM(SLSimpleurdf, optionsSimple);
nqS = SLSimple.getNumPositions();
nvS = SLSimple.getNumVelocities();


%% Define NLP

qS_inds = reshape(1:(N*nqS), nqS, N);
vS_inds = reshape(N*nqS + (1:(N*nvS)), nvS, N);

num_vars = (nqS + nvS)*N;
x_name = cell(num_vars, 1);

for i = 1:N
    for j = 1:nqS
        x_name{(i-1)*nqS+j} = sprintf('qS%d[%d]',j,i);
    end
    for j = 1:nvS
        x_name{N*nqS+(i-1)*nvS+j} = sprintf('vS%d[%d]',j,i);
    end
end

nlp = NonlinearProgram(num_vars, x_name);

nlp = nlp.setSolver('snopt');
nlp = nlp.setSolverOptions('snopt','MajorIterationsLimit',10000);
nlp = nlp.setSolverOptions('snopt','MinorItxerationsLimit',200000);
nlp = nlp.setSolverOptions('snopt','IterationsLimit',5000000);
nlp = nlp.setSolverOptions('snopt','SuperbasicsLimit',1000);

tol = 1e-6;
nlp = nlp.setSolverOptions('snopt','MajorOptimalityTolerance',tol);
nlp = nlp.setSolverOptions('snopt','MinorOptimalityTolerance',tol);
nlp = nlp.setSolverOptions('snopt','MajorFeasibilityTolerance',tol);
nlp = nlp.setSolverOptions('snopt','MinorFeasibilityTolerance',tol);
nlp = nlp.setSolverOptions('snopt','constraint_err_tol',tol);


%% Add objective

for i = 1:N   
    nlp = nlp.addCost(FunctionHandleObjective(nqS+nvS, @(qS, vS)objective_fun(...
        qS, vS, xyz_foot(:,i), vxyz_foot(:,i))), {qS_inds(:, i); vS_inds(:,i)});    
end

%% Initialize and Solve

qS0 = bsxfun(@plus, XI(SFB_INDS, :), SLSimple.q0); 
vS0 = XI(nq+SFB_INDS, :); 
xS0 = [qS0(:); vS0(:)]; 
[xf, objval] = nlp.solve(xS0);


%% Objective
    function [f, df] = objective_fun(qS, vS, xyz_foot, vxyz_foot)
        
        xin = [qS; vS];
        [f,df] = objective(xin, xyz_foot, vxyz_foot);
%         fprintf('Objective: %f \r', max(f));
        
%         step = 1e-6;
%         df_fd = zeros(size(df));
%         dxin = step*eye(length(xin));
%         for k = 1:length(xin)
%             df_fd(:,k) = (objective(xin+dxin(:,k), xyz_foot, vxyz_foot) - ...
%                 objective(xin-dxin(:,k), xyz_foot, vxyz_foot))/(2*step);
%         end
%         
%         disp('Objective Derivative Error:');
%         disp(max(abs(df_fd(:)-df(:))));
        
    end


    function [f, df] = objective(xin, xyz_foot, vxyz_foot)
        
        qS = xin(1:nqS);
        vS = xin(nqS+(1:nvS));
        dim = numel(xyz_foot);          

        kinsol = SLSimple.doKinematics(qS, vS, struct('compute_gradients', true));
        [xyzS_foot, JSfoot, dJSfoot] = SLSimple.forwardKin(kinsol, SLSimple.findLinkId('L2'), ...
            [0 0 -14.988382167532292]');
        vxyzS_foot = JSfoot*vS;         
        
        x = [xyzS_foot; vxyzS_foot]; 
        xd = [xyz_foot; vxyz_foot];
        f = (1/2)*(x - xd)'*(x - xd);
        df = (x - xd)'*[JSfoot, zeros(dim, nvS);  ...
            reshape(reshape(dJSfoot, dim*nqS, nqS)*vS, dim, nqS), JSfoot];%
    end
end

 