function FlyStraight()

p = McFoamy();

x0 = [0 0 1]';
xf = [1 0 1]';

N = 2;
prog = DircolTrajectoryOptimization(p,N,[0 1]);
prog = prog.addStateConstraint(ConstantConstraint([x0; 10; 0; 0; 0; 0; 0]), 1, [1 2 3 8 9 10 11 12 13]);
prog = prog.addStateConstraint(QuadraticConstraint(1,1,2*eye(4),zeros(4,1)), 1:N, [4 5 6 7]);
prog = prog.addInputConstraint(BoundingBoxConstraint([0 -.5 -.5 -.5]', [700 .5 .5 .5]), 1:N);
prog = prog.addRunningCost(@cost);

traj_init.x = PPTrajectory(foh([0 .1], [[x0; 0; 1; 0; 0; 5; 0; 0; 0; 0; 0], [xf; 0; 1; 0; 0; 5; 0; 0; 0; 0; 0]]));
traj_init.u = PPTrajectory(foh([0 .1], [[500; 0; 0; 0], [500; 0; 0; 0]]));

[xtraj,utraj,~,~,info] = solveTraj(prog,.1,traj_init);

    function [g,dg] = cost(dt,x,u)
        R = diag([1 0 0 0]);
        Q = .1*diag([0 0 0 0 0 0 0 0 1 1 1 1 1]);
        g = x'*Q*x + u'*R*u;
        dg = [0, 2*x'*Q, 2*u'*R];
    end

end

