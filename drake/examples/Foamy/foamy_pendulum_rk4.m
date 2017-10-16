function [xn, xdot] = foamy_pendulum_rk4(x,u,dt)

xdot1 = foamy_pendulum_dynamics_mex(0,x,u);
if x(3) < 0.01
    xdot1 = clip(xdot1);
end

xm1 = x + (0.5*dt)*xdot1;
xm1 = normalize(xm1);

xdot2 = foamy_pendulum_dynamics_mex(0,xm1,u);
if xm1(3) < 0.01
    xdot2 = clip(xdot2);
end
xm2 = x + (0.5*dt)*xdot2;
xm2 = normalize(xm2);

xdot3 = foamy_pendulum_dynamics_mex(0,xm2,u);
if xm2(3) < 0.01
    xdot3 = clip(xdot3);
end

xm3 = x + dt*xdot3;
xm3 = normalize(xm3);

xdot4 = foamy_pendulum_dynamics_mex(0,xm3,u);
if xm3(3) < 0.01
    xdot4 = clip(xdot4);
end

xdot = (1/6)*(xdot1 + 2*xdot2 + 2*xdot3 + xdot4);
xdot = clip(xdot);
xn = x + dt*xdot;

xn = normalize(xn);

end

function xdotc = clip(xdot)
    %Make sure we don't go through the ground
    xdotc = xdot;
    xdotc(3) = max(0,xdotc(3));
    xdotc(10) = xdotc(3);
    xdotc(17) = max(0,xdotc(17));
    xdotc(23) = xdotc(17);
    xdotc(11) = 0;
end

function xn = normalize(x)
    xn = x;
    xn(4:7) = xn(4:7)/sqrt(xn(4:7)'*xn(4:7));
    xn(11:14) = xn(11:14)/sqrt(xn(11:14)'*xn(11:14));
end

