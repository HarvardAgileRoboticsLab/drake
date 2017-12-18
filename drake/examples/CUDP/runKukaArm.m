function runKukaArm2()
    % get plant and visualizer
    options.with_weight = true;
    options.with_shelf_and_boxes = true;
    p = KukaArm(options);
    v = p.constructVisualizer;
    v.display_dt = .05;
    
    % get constrants
    nx = p.getNumStates;
    nq = p.getNumPositions;
    nu = p.getNumInputs;
    tf0 = 4;
    N = tf0*100;
    tourqueLimit = 200;
    SigmaScale = 1e-4;
    tols = [[1e-2,5e-1];[1e-3,1e-2];[5e-5,5e-3]];
    Q = eye(nq*2);
    R = .0001*eye(nu);
    Qf = 1000*eye(nq*2);

    % create initial and final states and default trajectory
    q0 = [0;-0.683;0;1.77;0;0.88;-1.57];
    x0 = [q0;zeros(nq,1)];
    xf = [0;0;0;-pi/4;0;pi/4;pi/2;zeros(7,1)];
    traj_init.x = PPTrajectory(foh([0,tf0],[x0,xf]));
    traj_init.u = ConstantTrajectory(zeros(nu,N-1));
    
    % prep the visualizer
    v.draw(0,x0);

    % set up initial optimization problem
    options.debug = true;
    options.tolFun = 10;
    options.tolFunExit = 1e-2;
    options.SigmaScale = SigmaScale;
    options.muMax = 1e30;
    options.phiFactor = 10;
    options.muFactor = 10;
    options.ILambdaDecay = 1;
    % penalty iLQR
    options.penaltyOnly = true;
    options.tolConstr = tols(1,1);
    options.tolLambdaUpdate = tols(1,2);
    options.derivative_method = 1;
    traj_opt1a = setUpTraj(p,N,tf0,xf,nx,nq,nu,tourqueLimit,Q,R,Qf,options);
    options.tolConstr = tols(2,1);
    options.tolLambdaUpdate = tols(2,2);
    traj_opt1b = setUpTraj(p,N,tf0,xf,nx,nq,nu,tourqueLimit,Q,R,Qf,options);
    options.tolConstr = tols(3,1);
    options.tolLambdaUpdate = tols(3,2);
    traj_opt1c = setUpTraj(p,N,tf0,xf,nx,nq,nu,tourqueLimit,Q,R,Qf,options);

    % penalty UDP
    options.tolConstr = tols(1,1);
    options.tolLambdaUpdate = tols(1,2);
    options.derivative_method = 2;
    traj_opt2a = setUpTraj(p,N,tf0,xf,nx,nq,nu,tourqueLimit,Q,R,Qf,options);
    options.tolConstr = tols(2,1);
    options.tolLambdaUpdate = tols(2,2);
    traj_opt2b = setUpTraj(p,N,tf0,xf,nx,nq,nu,tourqueLimit,Q,R,Qf,options);
    options.tolConstr = tols(3,1);
    options.tolLambdaUpdate = tols(3,2);
    traj_opt2c = setUpTraj(p,N,tf0,xf,nx,nq,nu,tourqueLimit,Q,R,Qf,options);

    % AL iLQR
    options.penaltyOnly = false;
    options.derivative_method = 1;
    options.tolConstr = tols(1,1);
    options.tolLambdaUpdate = tols(1,2);
    traj_opt3a = setUpTraj(p,N,tf0,xf,nx,nq,nu,tourqueLimit,Q,R,Qf,options);
    options.tolConstr = tols(2,1);
    options.tolLambdaUpdate = tols(2,2);
    traj_opt3b = setUpTraj(p,N,tf0,xf,nx,nq,nu,tourqueLimit,Q,R,Qf,options);
    options.tolConstr = tols(3,1);
    options.tolLambdaUpdate = tols(3,2);
    traj_opt3c = setUpTraj(p,N,tf0,xf,nx,nq,nu,tourqueLimit,Q,R,Qf,options);

    % AL UDP
    options.derivative_method = 2;
    options.tolConstr = tols(1,1);
    options.tolLambdaUpdate = tols(1,2);
    traj_opt4a = setUpTraj(p,N,tf0,xf,nx,nq,nu,tourqueLimit,Q,R,Qf,options);
    options.tolConstr = tols(2,1);
    options.tolLambdaUpdate = tols(2,2);
    traj_opt4b = setUpTraj(p,N,tf0,xf,nx,nq,nu,tourqueLimit,Q,R,Qf,options);
    options.tolConstr = tols(3,1);
    options.tolLambdaUpdate = tols(3,2);
    traj_opt4c = setUpTraj(p,N,tf0,xf,nx,nq,nu,tourqueLimit,Q,R,Qf,options);

    % Solve
    % penalty iLQR
    tic;
%     [x1a,u1a,d1a] = traj_opt1a.solveTraj(tf0,traj_init);
    t1a = toc;
    tic;
%     [x1b,u1b,d1b] = traj_opt1b.solveTraj(tf0,traj_init);
    t1b = toc;
    tic;
%     [x1c,u1c,d1c] = traj_opt1c.solveTraj(tf0,traj_init);
    t1c = toc;
    % penalty UDP
    tic;
    [x2a,u2a,d2a] = traj_opt2a.solveTraj(tf0,traj_init);
    t2a = toc;
    tic;
%     [x2b,u2b,d2b] = traj_opt2b.solveTraj(tf0,traj_init);
    t2b = toc;
    tic;
%     [x2c,u2c,d2c] = traj_opt2c.solveTraj(tf0,traj_init);
    t2c = toc;
    % AL iLQR
    tic;
%     [x3a,u3a,d3a] = traj_opt3a.solveTraj(tf0,traj_init);
    t3a = toc;
    tic;
%     [x3b,u3b,d3b] = traj_opt3b.solveTraj(tf0,traj_init);
    t3b = toc;
    tic;
%     [x3c,u3c,d3c] = traj_opt3c.solveTraj(tf0,traj_init);
    t3c = toc;
%     % AL UDP
    tic;
%     [x4a,u4a,d4a] = traj_opt4a.solveTraj(tf0,traj_init);
    t4a = toc;
    tic;
%     [x4b,u4b,d4b] = traj_opt4b.solveTraj(tf0,traj_init);
    t4b = toc;
    tic;
%     [x4c,u4c,d4c] = traj_opt4c.solveTraj(tf0,traj_init);
    t4c = toc;
    
    keyboard;
    v.playback(x4c);

    %--------- Cost + Constraint Functions ---------%
    function [g,dg,ddg] = cost(dt,x,u,xf,nq,nu,Q,R);
        g = (x-xf)'*Q*(x-xf) + u'*R*u;
        dg = [0, 2*(x-xf)'*Q,2*u'*R];
        ddg = blkdiag(0,2*Q,2*R);
    end

    function [h,dh,ddh] = finalCost(tf,x,xf,nq,Qf)
        h = (x-xf)'*Qf*(x-xf);
        dh = [0, 2*(x-xf)'*Qf];
        ddh = blkdiag(0,2*Qf);
    end

    function [c,dc] = torqueLimit(k,u,tourqueLimit)
        c = [u - tourqueLimit; -u - tourqueLimit];
        dc = [eye(7); -eye(7)];
    end

    function [c,dc] = jointLimits(k,x,p)
        c = [x(1:7) - p.joint_limit_max; -x(1:7) + p.joint_limit_min];
        dc = [eye(7),zeros(7,7);-eye(7),zeros(7,7)];
    end

    function [c,dc] = xfConstraint(k,x,p,N,nx,nq)
        if k == N
            q = x(1:nq);
            kinsol = doKinematics(p, q);
            opt.rotation_type = 1;
            [pl,Jl] = forwardKin(p,kinsol,p.findLinkId('iiwa_link_ee'),[0;0;0],opt);

            c = pl-[0.4;0.0;1.03;pi/2;0;pi/2];
            dc = [Jl,zeros(length(c),nq)];
        else
            c = zeros(6,1);
            dc = zeros(6,nx);
        end
    end

    function [c,dc] = noContactConstraint(k,x,N,p,nq)
        if k < N
            q = x(1:nq);
            kinsol = doKinematics(p, q);
            [phi,~,~,~,~,~,~,~,J] = p.contactConstraints(kinsol);
            c = -1*phi;
            dc = -1*[J,zeros(length(c),nq)];
        else
            c = zeros(4,1);
            dc = zeros(4,2*nq);
        end
    end

    function [c,dc] = jointLimitsPlusContactConstraint(k,x,N,p,nq)
        [c1,dc1] = jointLimits(k,x,p);
        [c2,dc2] = noContactConstraint(k,x,N,p,nq);
        c = [c1;c2];
        dc = [dc1;dc2];
    end

    %--------------- Plan Visualizer Hook ---------------%
    function traj_opt = addPlanVisualizer(p,traj_opt,N)
      % spew out an lcmgl visualization of the trajectory.  intended to be
      % used as a callback (fake objective) in the direct trajectory
      % optimization classes

      if ~checkDependency('lcmgl')
        warning('lcmgl dependency is missing.  skipping visualization'); 
        return;
      end
      lcmgl = drake.util.BotLCMGLClient(lcm.lcm.LCM.getSingleton(), 'KukaArmPlan');
      
      traj_opt = traj_opt.addPlanVisualizer(@(x)visualizePlan(x,lcmgl,p,N),[1,2,3,4,5,6,7,8,9,10,11,12,13,14]);
      
      function visualizePlan(x,lcmgl,p,N)
        handx = zeros(1,N);
        handy = zeros(1,N);
        handz = zeros(1,N);
        for k = 1:N
            hand_pos = handPos(p,x(:,k));
            handx(:,k) = hand_pos(1);
            handy(:,k) = hand_pos(2);
            handz(:,k) = hand_pos(3);
        end
        lcmgl.glColor3f(1, 0, 0);
        lcmgl.glPointSize(3);
        lcmgl.points(handx,handy,handz);
        lcmgl.glColor3f(.5, .5, 1);
        lcmgl.plot3(handx,handy,handz);
        lcmgl.switchBuffers;
      end
      
      function [hand_pos, dHand_pos, ddHand_pos] = handPos(p,x)
        if nargout > 1
            q = x(1:p.num_positions);
            qd = x(p.num_positions+1:p.num_positions+p.num_velocities);
            kinsol_options.compute_gradients = true;
            kinsol = p.doKinematics(q,qd,kinsol_options);
            hand_body = p.findLinkId('iiwa_link_7');
            pos_on_hand_body = [0;0;0];
            [hand_pos, dHand_pos, ddHand_pos] = p.forwardKin(kinsol,hand_body,pos_on_hand_body);
        else
            q = x(1:p.num_positions);
            qd = x(p.num_positions+1:p.num_positions+p.num_velocities);
            kinsol_options.compute_gradients = false;
            kinsol = p.doKinematics(q,qd,kinsol_options);
            hand_body = p.findLinkId('iiwa_link_7');
            pos_on_hand_body = [0;0;0.25];
            hand_pos = p.forwardKin(kinsol,hand_body,pos_on_hand_body);
        end
      end
    end

    %-------------------HELPERS---------------------------%
    function [traj_opt] = setUpTraj(p,N,tf0,xf,nx,nq,nu,tourqueLimit,Q,R,Qf,options)
        traj_opt = CUDPTrajectoryOptimization(p,N,[tf0 2*tf0],options);
        traj_opt = traj_opt.addRunningCost(@(dt,x,u)cost(dt,x,u,xf,nq,nu,Q,R));
        traj_opt = traj_opt.addFinalCost(@(tf,x)finalCost(tf,x,xf,nq,Qf));
        traj_opt = traj_opt.addInputConstraint(@(k,x)torqueLimit(k,x,tourqueLimit), 2*nu);
        traj_opt = traj_opt.addStateConstraint(@(k,x)xfConstraint(k,x,p,N,nx,nq), 6, 1);
        traj_opt = traj_opt.addStateConstraint(@(k,x)jointLimitsPlusContactConstraint(k,x,N,p,nq), 2*nq + 4, 0);
        traj_opt = addPlanVisualizer(p,traj_opt,N);
    end
end
