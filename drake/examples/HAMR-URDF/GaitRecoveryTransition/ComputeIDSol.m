function [xtraj, utraj, vtraj, z, F, info] = ComputeIDSol(hamr, actuators, trajd, v, params)

% unpack params
NSAMP = params.NSAMP;
T = params.T;
h = T/NSAMP;

% dimensions
nq = params.nq;
nv = params.nv;
nu = params.nu;
nl = params.nl;
qa = params.qa;
nqa = params.nqa;

% build trajectory optimization
% options.periodic = true; 
x0 = [zeros(nq,1); zeros(nv,1)];
traj_opt = DirtranTrajectoryOptimization(hamr, NSAMP, T); %, options);
traj_init.x = PPTrajectory(foh([0,T],[x0,x0]));

% initial state
% traj_opt = traj_opt.addStateConstraint(ConstantConstraint(trajd(1, 2:end)'), 1, qa);

% add periodic state constraint
QXperiodic = [-eye(nq+nv), eye(nq+nv)];       
periodic_state_cnst = LinearConstraint(zeros(nq+nv,1), zeros(nq+nv,1),QXperiodic);
periodic_state_cnst = periodic_state_cnst.setName('periodicity_state');
traj_opt = traj_opt.addStateConstraint(periodic_state_cnst,{[1, NSAMP]});
 
% add periodic input constraint
QUperiodic = (1/2)*[-eye(nu), -eye(nu), eye(nu), eye(nu)];%blkdiag([-eye(nu), zeros(nu), eye(nu), zeros(nu)], [-eye(nu), eye(nu), eye(nu), -eye(nu)]);
periodic_input_cnst = LinearConstraint(zeros(nu,1),zeros(nu,1),QUperiodic);
periodic_input_cnst = periodic_input_cnst.setName('periodicity_input');
traj_opt = traj_opt.addInputConstraint(periodic_input_cnst,{[1, 2, NSAMP-1, NSAMP]});


% add tracking constraint
td = linspace(0, T, size(trajd, 1));
qqd = interp1(td, trajd(:, 2:end), linspace(0, T, NSAMP)); 
% figure(1); plot(qqd(:,1), '*-')
% vvd = diff(qqd)/h; 
% vvd = [vvd; vvd(end, :)];
% track nominal trajectory (loose)
for ind = 2:(NSAMP-1)
    traj_opt = traj_opt.addCost(FunctionHandleObjective(nqa,@(x)tracking_cost(x, qqd(ind,:)'),1), ...
        traj_opt.x_inds(qa,ind));
end

% try to get smooth u
for ind = 1:NSAMP-1
    traj_opt = traj_opt.addCost(FunctionHandleObjective(2*(nu-nl),@(u1, u2)udiff_cost(u1, u2)), ...
        {traj_opt.u_inds(1:nu-nl,ind); traj_opt.u_inds(1:nu-nl,ind+1)});
end


% add loop constraints
loop_constraint = FunctionHandleConstraint(zeros(nl,1), zeros(nl,1), nq+nv, @(x)loop_const(hamr, x));
traj_opt = traj_opt.addStateConstraint(loop_constraint, 1:(NSAMP));

%add quadratic input cost
% traj_opt = traj_opt.addRunningCost(@running_cost);

% add display function
traj_opt = traj_opt.addTrajectoryDisplayFunction(@displayTraj);


% ----- Solver options -----
traj_opt = traj_opt.setSolver('snopt');
traj_opt = traj_opt.setSolverOptions('snopt','MajorIterationsLimit',10000);
traj_opt = traj_opt.setSolverOptions('snopt','MinorIterationsLimit',200000);
traj_opt = traj_opt.setSolverOptions('snopt','IterationsLimit',10000000);
traj_opt = traj_opt.setSolverOptions('snopt','SuperbasicsLimit',1000);

traj_opt = traj_opt.setSolverOptions('snopt','MajorOptimalityTolerance',1e-5);
traj_opt = traj_opt.setSolverOptions('snopt','MinorOptimalityTolerance',1e-5);
traj_opt = traj_opt.setSolverOptions('snopt','MajorFeasibilityTolerance',1e-5);
traj_opt = traj_opt.setSolverOptions('snopt','MinorFeasibilityTolerance',1e-5);
traj_opt = traj_opt.setSolverOptions('snopt','constraint_err_tol',1e-5);

% solve
[xtraj,utraj,z,F,info] = traj_opt.solveTraj(T,traj_init);

% compute actuator voltages
ttopt = xtraj.getBreaks();
xxopt = xtraj.eval(ttopt);
uuopt = utraj.eval(ttopt);

vvopt = zeros(nu-nl, numel(ttopt));
for j = 1:numel(ttopt)
    vvj = ipzt_fun(hamr, actuators, ttopt(j), xxopt(qa,j), uuopt(1:(nu-nl),j));
    vvopt(:,j) = vvj(1:(nu-nl));
end

vvopt(vvopt < actuators.dummy_bender(1).dp.Vg, j) = actuators.dummy_bender(1).dp.Vg;
vvopt(vvopt > actuators.dummy_bender(1).dp.Vb, j) = actuators.dummy_bender(1).dp.Vb;
vtraj = PPTrajectory(foh(ttopt, vvopt));



% ---------------- Helper Functions ---------------------------------

    function [c, dc] = tracking_cost(x, xd)
        Q = eye(nqa);
        c = (1/2)*(x-xd)'*Q*(x-xd);
        dc = (x-xd)'*Q;
        %         fprintf('Tracking Cost: %f \r', c)
    end

%     function [c, dc] = running_cost(dt, x, u)
%         R = eye(nu);
%         u0 = (u-actuators.dummy_bender(1).dp.Vb/2); 
%         c = (1/2)*u0'*R*u0;
%         dc = [0, zeros(1, nq+nv), u0'*R];
%     end

    function [f, df] = udiff_cost(u1, u2)
        udiff = u1-u2;
        a = 1; 
        f = (1/2)*(udiff'*udiff);
        I = eye(length(u1));
        df = a*[udiff'*I,-udiff'*I];
    end

    function [f, df] = loop_const(obj, x)
        q = x(1:nq);
%         qd = x(nq+(1:nv)); 
        [k, K, ~] = obj.positionConstraints(q);
%         dK = reshape(dK(obj.valid_loops, :)', nq, nl*nq)';

        f = k(obj.valid_loops); % K(obj.valid_loops,:)*qd]; 
        df = [K(obj.valid_loops, :), zeros(nl, nv)]; ; ...
            %kron(qd', eye(nl))*dK, K(obj.valid_loops,:)];
    end

    function displayTraj(h,x,u)
%         disp('Displaying Trajectory...')
        h = h/1e3;
        ts = [0;cumsum(h)];
        for k=1:length(ts)
            v.drawWrapper(0,x(:,k));
            pause(h(1));
        end
    end

end