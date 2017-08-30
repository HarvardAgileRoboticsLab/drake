function [xn, xdot] = foamy_rk4(x,u,dt)

xdot1 = foamy_dynamics_mex(0,x,u);
if x(3) < 0.01
    xdot1 = clip(xdot1);
end

xm1 = x + (0.5*dt)*xdot1;
xm1(4:7) = xm1(4:7)/sqrt(xm1(4:7)'*xm1(4:7));

xdot2 = foamy_dynamics_mex(0,xm1,u);
if xm1(3) < 0.01
    xdot2 = clip(xdot2);
end
xm2 = x + (0.5*dt)*xdot2;
xm2(4:7) = xm2(4:7)/sqrt(xm2(4:7)'*xm2(4:7));

xdot3 = foamy_dynamics_mex(0,xm2,u);
if xm2(3) < 0.01
    xdot3 = clip(xdot3);
end

xm3 = x + dt*xdot3;
xm3(4:7) = xm3(4:7)/sqrt(xm3(4:7)'*xm3(4:7));

xdot4 = foamy_dynamics_mex(0,xm3,u);
if xm3(3) < 0.01
    xdot4 = clip(xdot4);
end

xdot = (1/6)*(xdot1 + 2*xdot2 + 2*xdot3 + xdot4);
xn = x + dt*xdot;

xn(3) = max(0,xn(3));
if xn(3) == 0
    xdot = clip(xdot);
    xn(10) = max(0,xn(10));
    xn(11) = 0;
end

end

function xdotc = clip(xdot)
    %Make sure we don't go through the ground
    xdotc = xdot;
    xdotc(3) = max(0,xdotc(3));
    xdotc(10) = max(0,xdotc(10));
    xdotc(11) = 0;
end

