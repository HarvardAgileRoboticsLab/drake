function [xn, xdot] = foamy_rk4(x,u,dt)

xdot1 = foamy_dynamics_mex(0,x,u);
[~,xdot1] = clip(x,xdot1);
xm1 = x + (0.5*dt)*xdot1;
xm1 = clip(xm1);

xdot2 = foamy_dynamics_mex(0,xm1,u);
[~,xdot2] = clip(xm1,xdot2);
xm2 = x + (0.5*dt)*xdot2;
xm2 = clip(xm2);

xdot3 = foamy_dynamics_mex(0,xm2,u);
[~,xdot3] = clip(xm2,xdot3);
xm3 = x + dt*xdot3;

xdot4 = foamy_dynamics_mex(0,xm3,u);
[~,xdot4] = clip(xm3,xdot4);

xdot = (1/6)*(xdot1 + 2*xdot2 + 2*xdot3 + xdot4);
xn = x + dt*xdot;
[xn,xdot] = clip(xn,xdot);

end

function [xc,xdotc] = clip(x,varargin)

if nargin == 2
    xdot = varargin{1};
else
    xdot = zeros(13,1);
end

xc = x;
xdotc = xdot;

%Make sure quaternion is normalized
xc(4:7) = x(4:7)/sqrt(x(4:7)'*x(4:7));

%Make sure we don't go through the ground
xc(3) = max(0,x(3));
if(xc(3) == 0)
    xc(10) = max(0,x(10));
    xdotc(3) = max(0,xdot(3));
    xdotc(10) = max(0,xdot(10));
end
end
