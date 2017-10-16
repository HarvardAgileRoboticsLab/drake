function [xn, xdot] = foamy_pendulum_midpoint(x,u,dt)

xdot = foamy_pendulum_dynamics_mex(0,x,u);

%Make sure we don't go through the ground
if(x(3) < 0.05)
    xdot(3) = max(0,xdot(3));
    xdot(10) = xdot(3);
    xdot(17) = max(0,xdot(17));
    xdot(23) = xdot(17);
    xdot(18) = 0;
    xdot(24) = 0;
end

xm = x + (0.5*dt)*xdot;

%Make sure quaternion is normalized
xm(4:7) = xm(4:7)/sqrt(xm(4:7)'*xm(4:7));
xm(11:14) = xm(11:14)/sqrt(xm(11:14)'*xm(11:14));

xdot = foamy_pendulum_dynamics_mex(0,xm,u);

%Make sure we don't go through the ground
if(x(3) < 0.05)
    xdot(3) = max(0,xdot(3));
    xdot(10) = xdot(3);
    xdot(17) = max(0,xdot(17));
    xdot(23) = xdot(17);
    xdot(18) = 0;
    xdot(24) = 0;
end

xn = x + dt*xdot;

%Make sure quaternion is normalized
xn(4:7) = xn(4:7)/sqrt(xn(4:7)'*xn(4:7));
xn(11:14) = xn(11:14)/sqrt(xn(11:14)'*xn(11:14));


end

