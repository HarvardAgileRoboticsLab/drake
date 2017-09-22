function [xn, xdot] = foamy_midpoint(x,u,dt)

xdot = foamy_dynamics_mex(0,x,u);

%Make sure we don't go through the ground
if(x(3) < 0.05)
    %R = qtoR(x(4:7));
    %P = R*diag([1 0 1])*R';
    %xdot(1:3) = P*xdot(1:3);
    xdot(3) = max(0,xdot(3));
    %xdot(8:10) = P*xdot(8:10);
    xdot(10) = max(0,xdot(10));
    xdot(11) = 0;
end

xm = x + (0.5*dt)*xdot;

%Make sure quaternion is normalized
xm(4:7) = xm(4:7)/sqrt(xm(4:7)'*xm(4:7));

xdot = foamy_dynamics_mex(0,xm,u);

%Make sure we don't go through the ground
if(x(3) < 0.05)
    %R = qtoR(xm(4:7));
    %P = R*diag([1 0 1])*R';
    %xdot(1:3) = P*xdot(1:3);
    xdot(3) = max(0,xdot(3));
    %xdot(8:10) = P*xdot(8:10);
    xdot(10) = max(0,xdot(10));
    xdot(11) = 0;
end

xn = x + dt*xdot;

%Make sure quaternion is normalized
xn(4:7) = xn(4:7)/sqrt(xn(4:7)'*xn(4:7));

%Make sure we don't go through the ground
xn(3) = max(0,xn(3));
if xn(3) < 0.05
    %R = qtoR(xm(4:7));
    %P = R*diag([1 0 1])*R';
    %xn(8:10) = P*xn(8:10);
    xn(10) = max(0,xn(10));
    xn(11) = 0;
end

    function R = qtoR(q)
        R = eye(3) + 2*hat(q(2:4))*(hat(q(2:4)) + q(1)*eye(3));
    end

    function C = hat(v)
        C = [0   -v(3) v(2); v(3)  0   -v(1); -v(2) v(1)  0];
    end

end

