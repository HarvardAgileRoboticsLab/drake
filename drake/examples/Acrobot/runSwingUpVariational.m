function [] = runSwingUpVariational()

options.ignore_self_collisions = true;
options.use_bullet = false;
p = PlanarRigidBodyManipulator('Acrobot.urdf',options);
v = p.constructVisualizer();

x0 = zeros(4,1);
xf = [pi 0 0 0]';
tf0 = 3;

N = 31;
prog = VariationalTrajectoryOptimization(p,N,[2 6]);
%prog = prog.addStateConstraint(ConstantConstraint(x0),1);
%prog = prog.addStateConstraint(ConstantConstraint(xf),N);
prog = prog.addPositionConstraint(ConstantConstraint(x0(1:2)),1);
prog = prog.addPositionConstraint(ConstantConstraint(x0(1:2)),2);
%prog = prog.addVelocityConstraint(ConstantConstraint(x0(3:4)),1);
%prog = prog.addInputConstraint(ConstantConstraint(0),1:N-1);
prog = prog.addPositionConstraint(ConstantConstraint(xf(1:2)),N-1);
prog = prog.addPositionConstraint(ConstantConstraint(xf(1:2)),N);

prog = prog.addRunningCost(@cost);
prog = prog.addFinalCost(@finalCost);

traj_init.x = PPTrajectory(foh([0,tf0],[x0,xf]));
traj_init.u = PPTrajectory(zoh(linspace(0,tf0,N), .1*randn(1,N)));

tic
[xtraj,utraj,z,F,info] = prog.solveTraj(tf0,traj_init);
toc

v.playback(xtraj, struct('slider', true));

    function [g,dg] = cost(h,x,u)
        R = 1;
        Q = diag([10 10 1 1]);
        xerr = x-xf;
        xerr(1,:) = mod(xerr(1,:)+pi,2*pi)-pi;
        g = h*(xerr'*Q*xerr + u'*R*u);
        dg = [g/h, 2*h*xerr'*Q, 2*h*u'*R];
    end

    function [h,dh] = finalCost(t,x)
        h = t;
        dh = [1,zeros(1,size(x,1))];
    end
end
