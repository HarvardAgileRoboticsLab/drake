function runQuad()

    % add all of the paths we need
    addpath('../Quadrotor');

    p = Quadrotor();
    % add the trees
    quadRadius = 0.5;
    treeSize = [0.5,0.5,2];
    treeLocations = [[0;1],[2;0],[2;1],[4;1],[6;1],[0.75;3],[3;3],[6;3],[7;3],[6;2],[6;4],[6;5],[3;4],[3;5],[4;5],[0;5],[1;7],[5;7],[7;7],[0;9],[2;9],[4;9],[6.5;9],[7;8],[3;8]];
    NTrees = max(size(treeLocations));
    for k = 1:NTrees
        p = addTree(p, treeSize, treeLocations(:,k), 0);
    end
    % construct the visualizer
    v = constructVisualizer(p);
    
    % set the initial and final conditions
    x0 = [0;0;0.5;zeros(9,1)];
    xf = [7;10;0.5;zeros(9,1)];
    
    % start the visualizer
    v.draw(0,double(x0));
    
    % prep for the solver
    tf0 = 4;
    N = 120;
    tourqueLimit = 10;
    traj_init.x = PPTrajectory(foh([0,tf0],[double(x0),double(xf)]));
    
    Qf = 1000*eye(12);
    Q = .1*eye(12);
    R = .01*eye(4);
    SigmaScale = 1e-4;
    tols = [[1e-2,1e-1];[1e-4,1e-2];[1e-6,5e-3]];
   
    % set up solvers
    options.debug = true;
    options.SigmaScale = SigmaScale;
    options.muMax = 1e30;
    options.phiFactor = 2;
    options.muFactor = 10;
    options.ILambdaDecay = 1;
    % penalty iLQR
    options.penaltyOnly = true;
    options.derivative_method = 1;
    options.tolConstr = tols(1,1);
    options.tolLambdaUpdate = tols(1,2);
    traj_opt1a = setUpTrajOpt(p,N,options,xf,tourqueLimit,Qf,Q,R,NTrees,quadRadius,treeSize,treeLocations);
    options.tolConstr = tols(2,1);
    options.tolLambdaUpdate = tols(2,2);
    traj_opt1b = setUpTrajOpt(p,N,options,xf,tourqueLimit,Qf,Q,R,NTrees,quadRadius,treeSize,treeLocations);
    options.tolConstr = tols(3,1);
    options.tolLambdaUpdate = tols(3,2);
    traj_opt1c = setUpTrajOpt(p,N,options,xf,tourqueLimit,Qf,Q,R,NTrees,quadRadius,treeSize,treeLocations);

    % penalty UDP
    options.tolConstr = tols(1,1);
    options.tolLambdaUpdate = tols(1,2);
    options.derivative_method = 2;
    traj_opt2a = setUpTrajOpt(p,N,options,xf,tourqueLimit,Qf,Q,R,NTrees,quadRadius,treeSize,treeLocations);
    options.tolConstr = tols(2,1);
    options.tolLambdaUpdate = tols(2,2);
    traj_opt2b = setUpTrajOpt(p,N,options,xf,tourqueLimit,Qf,Q,R,NTrees,quadRadius,treeSize,treeLocations);
    options.tolConstr = tols(3,1);
    options.tolLambdaUpdate = tols(3,2);
    traj_opt2c = setUpTrajOpt(p,N,options,xf,tourqueLimit,Qf,Q,R,NTrees,quadRadius,treeSize,treeLocations);

    % AL iLQR
    options.penaltyOnly = false;
    options.derivative_method = 1;
    options.tolConstr = tols(1,1);
    options.tolLambdaUpdate = tols(1,2);
    traj_opt3a = setUpTrajOpt(p,N,options,xf,tourqueLimit,Qf,Q,R,NTrees,quadRadius,treeSize,treeLocations);
    options.tolConstr = tols(2,1);
    options.tolLambdaUpdate = tols(2,2);
    traj_opt3b = setUpTrajOpt(p,N,options,xf,tourqueLimit,Qf,Q,R,NTrees,quadRadius,treeSize,treeLocations);
    options.tolConstr = tols(3,1);
    options.tolLambdaUpdate = tols(3,2);
    traj_opt3c = setUpTrajOpt(p,N,options,xf,tourqueLimit,Qf,Q,R,NTrees,quadRadius,treeSize,treeLocations);

    % AL UDP
    options.derivative_method = 2;
    options.tolConstr = tols(1,1);
    options.tolLambdaUpdate = tols(1,2);
    traj_opt4a = setUpTrajOpt(p,N,options,xf,tourqueLimit,Qf,Q,R,NTrees,quadRadius,treeSize,treeLocations);
    options.tolConstr = tols(2,1);
    options.tolLambdaUpdate = tols(2,2);
    traj_opt4b = setUpTrajOpt(p,N,options,xf,tourqueLimit,Qf,Q,R,NTrees,quadRadius,treeSize,treeLocations);
    options.tolConstr = tols(3,1);
    options.tolLambdaUpdate = tols(3,2);
    traj_opt4c = setUpTrajOpt(p,N,options,xf,tourqueLimit,Qf,Q,R,NTrees,quadRadius,treeSize,treeLocations);


    % solve
    % penalty iLQR
    tic;
    [x1a,u1a,d1a] = traj_opt1a.solveTraj(tf0,traj_init);
    t1a = toc;
    tic;
    [x1b,u1b,d1b] = traj_opt1b.solveTraj(tf0,traj_init);
    t1b = toc;
    tic;
    [x1c,u1c,d1c] = traj_opt1c.solveTraj(tf0,traj_init);
    t1c = toc;
    % penalty UDP
    tic;
    [x2a,u2a,d2a] = traj_opt2a.solveTraj(tf0,traj_init);
    t2a = toc;
    tic;
    [x2b,u2b,d2b] = traj_opt2b.solveTraj(tf0,traj_init);
    t2b = toc;
    tic;
    [x2c,u2c,d2c] = traj_opt2c.solveTraj(tf0,traj_init);
    t2c = toc;
    % AL iLQR
    tic;
    [x3a,u3a,d3a] = traj_opt3a.solveTraj(tf0,traj_init);
    t3a = toc;
    tic;
    [x3b,u3b,d3b] = traj_opt3b.solveTraj(tf0,traj_init);
    t3b = toc;
    tic;
    [x3c,u3c,d3c] = traj_opt3c.solveTraj(tf0,traj_init);
    t3c = toc;
    % AL UDP
    tic;
    [x4a,u4a,d4a] = traj_opt4a.solveTraj(tf0,traj_init);
    t4a = toc;
    tic;
    [x4b,u4b,d4b] = traj_opt4b.solveTraj(tf0,traj_init);
    t4b = toc;
    tic;
    [x4c,u4c,d4c] = traj_opt4c.solveTraj(tf0,traj_init);
    t4c = toc;
    
    keyboard;
    v.playback(x4c);
   
    function [g,dg,ddg] = cost(dt,x,u,xf,Q,R)
        g = (x-xf)'*Q*(x-xf) + u'*R*u;
        dg = [0, 2*(x-xf)'*Q,2*u'*R];
        ddg = blkdiag(0,2*Q,2*R);
    end

    function [h,dh,ddh] = finalCost(tf,x,xf,Qf)
        h = (x-xf)'*Qf*(x-xf);
        dh = [0, 2*(x-xf)'*Qf];
        ddh = blkdiag(0,2*Qf);
    end

    function [c,dc] = torqueLimit(k,u,tourqueLimit)
        c = [u - tourqueLimit; -u - tourqueLimit];
        dc = [eye(4); -eye(4)];
    end

    function [c,dc] = xConstrEq(k,x,xf,N)
        if k < N
            c = zeros(12,1);
            dc = zeros(12,12);
        else
            c = x-xf;
            dc = eye(12);
        end
    end

    function [c,dc] = treeCnstr(k,x,NTrees,quadRadius,treeSize,treeLocations)
        radius = treeSize(1);
        c = zeros(NTrees,1);
        dc = zeros(NTrees,12);
        for k = 1:NTrees
            treeLoc = treeLocations(:,k);
            xyDiff = x(1:2) - treeLoc;
            xyDist2 = xyDiff' * xyDiff;
            value = (radius + quadRadius)^2 - xyDist2;
            c(k,:) = value;
            dc(k,:) = [-2*xyDiff',zeros(1,10)];
        end
    end

    function traj_opt = addPlanVisualizer(obj,traj_opt)
      % spew out an lcmgl visualization of the trajectory.  intended to be
      % used as a callback (fake objective) in the direct trajectory
      % optimization classes

        if ~checkDependency('lcmgl')
            warning('lcmgl dependency is missing.  skipping visualization'); 
            return;
        end
        lcmgl = drake.util.BotLCMGLClient(lcm.lcm.LCM.getSingleton(), 'QuadrotorPlan');
        traj_opt = traj_opt.addPlanVisualizer(@(x)visualizePlan(x,lcmgl),[1,2,3]);
        function visualizePlan(x,lcmgl)
            lcmgl.glColor3f(1, 0, 0);
            lcmgl.glPointSize(3);
            lcmgl.points(x(1,:),x(2,:),x(3,:));
            lcmgl.glColor3f(.5, .5, 1);
            lcmgl.plot3(x(1,:),x(2,:),x(3,:));
            lcmgl.switchBuffers;
        end
    end

    function [traj_opt] = setUpTrajOpt(p,N,options,xf,tourqueLimit,Qf,Q,R,NTrees,quadRadius,treeSize,treeLocations)
        traj_opt = CUDPTrajectoryOptimization(p,N,[2 6],options);
        traj_opt = traj_opt.addInputConstraint(@(k,u)torqueLimit(k,u,tourqueLimit), 8);
        traj_opt = traj_opt.addStateConstraint(@(k,x)xConstrEq(k,x,xf,N), 12, 1);
        traj_opt = traj_opt.addStateConstraint(@(k,x)treeCnstr(k,x,NTrees,quadRadius,treeSize,treeLocations), NTrees, 0);
        traj_opt = traj_opt.addRunningCost(@(dt,x,u)cost(dt,x,u,xf,Q,R));
        traj_opt = traj_opt.addFinalCost(@(tf,x)finalCost(tf,x,xf,Qf));
%         traj_opt = addPlanVisualizer(p,traj_opt);
    end

end