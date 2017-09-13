function [utraj,xtraj]=variationalSwingUp()

p = PlanarRigidBodyManipulator('Pendulum.urdf');
v = p.constructVisualizer();

x0 = [0;0];
xf = [pi;0];
tf0 = 4;
N = 25;

traj_opt = VariationalTrajectoryOptimization(p,N,[2 6]);
traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(x0(1)),1);
traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(x0(1)),2);
traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(xf(1)),N-1);
traj_opt = traj_opt.addPositionConstraint(ConstantConstraint(xf(1)),N);
traj_opt = traj_opt.addRunningCost(@cost);
traj_opt = traj_opt.addFinalCost(@finalCost);

traj_init.x = PPTrajectory(foh([0 tf0], [double(x0) double(xf)]));

    function [g,dg] = cost(dt,x,u);
        R = 10;
        g = (R*u).*u;
        
        if (nargout>1)
            dg = [zeros(1,3),2*u'*R];
        end
    end

    function [h,dh] = finalCost(tf,x)
        h = tf;
        if (nargout>1)
            dh = [1, zeros(1,2)];
        end
    end

info=0;
while (info~=1)
    tic
    [xtraj,utraj,z,F,info] = traj_opt.solveTraj(tf0,traj_init);
    toc
end

Q = diag([10 10]);
R = 1;
tv = tvlqr(p,xtraj,utraj,Q,R,10*Q);

clsys = feedback(p,tv);
cltraj = simulate(clsys,utraj.tspan,[0;0]);

v.playback(cltraj);

end