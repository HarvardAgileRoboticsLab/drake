function [xn, xdot] = foamy_midpoint(x,u,dt)

xdot = foamy_dynamics_mex(0,x,u);

%Make sure we don't go through the ground
if(x(3) < 0.01)
    xdot(3) = max(0,xdot(3));
    xdot(10) = max(0,xdot(10));
    xdot(11) = 0;
end

xm = x + (0.5*dt)*xdot;

%Make sure quaternion is normalized
xm(4:7) = xm(4:7)/sqrt(xm(4:7)'*xm(4:7));

xdot = foamy_dynamics_mex(0,xm,u);

%Make sure we don't go through the ground
if(x(3) < 0.01)
    xdot(3) = max(0,xdot(3));
    xdot(10) = max(0,xdot(10));
    xdot(11) = 0;
end

xn = x + dt*xdot;

%Make sure quaternion is normalized
xn(4:7) = xn(4:7)/sqrt(xn(4:7)'*xn(4:7));

%Make sure we don't go through the ground
xn(3) = max(0,xn(3));
if xn(3) == 0
    xn(10) = max(0,xn(10));
    xn(11) = 0;
end

end

