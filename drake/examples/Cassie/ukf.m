function [xhat, Phat] = ukf(f,h,x,y,u,P,Q,R,ts)

Nx = length(x);
Ny = length(y);
Np = Nx+2; %Number of sigma points
w = 1/Np; %Sigma point weights

%Make sure quaternion is normalized
x(4:7) = x(4:7)/norm(x(4:7));

%Generate sigma points
S = simplex_points(Nx);
L = chol(P,'lower');
S = L*S;
chi = zeros(Nx+1,Np+1);
chi(:,1) = x;
for j = 2:Np
    chi(:,j) = [x(1:3) + S(1:3,j); qmultiply(gtoq(S(4:6,j)),x(4:7)); x(8:13) + S(7:12,j)];
end

%Propagate sigma points forward one time step
x_p = zeros(Nx+1,Np);
for j = 1:Np
    x_p(:,j) = rkstep(chi(:,j), u, ts);
end

%Calculate predicted state
xbar = zeros(Nx+1,1);
M = zeros(4);
for j = 1:Np
    xbar(1:3) = xbar(1:3) + w*x_p(1:3,j);
    M = M + w*x_p(4:7,j)*x_p(4:7,j)';
    xbar(8:13) = xbar(8:13) + w*x_p(8:13,j);
end
%Quaternion average (from Markley et.al. 2007)
M = 4*M - eye(4);
qbar = x(4:7); e = 1; %initial guess
for l = 1:2 %Rayleigh quotient iteration
    qbar = (M - e*eye(4))\qbar;
    qbar = qbar/sqrt(qbar'*qbar);
    e = qbar'*M*qbar;
end
xbar(4:7) = qbar;

%Calculate predicted measurements
y_p = zeros(Ny,Np);
for j = 1:Np
    xdot_p = f(x_p(:,j),u);
    y_p(:,j) = h(x_p(:,j),xdot_p);
end
ybar = zeros(Ny,1);
for j = 1:Np
    ybar = ybar + w*y_p(:,j);
end

%Calculate predicted covariances
Pbar = Q;
Pyy = zeros(Ny,Ny);
Pxy = zeros(Nx,Ny);
for k = 1:Np
    rx = [x_p(1:3,k)-xbar(1:3); qtog(qmultiply(x_p(4:7,k),qconj(qbar))); x_p(8:13,k)-xbar(8:13)];
    ry = y_p(:,k)-ybar;
    Pbar = Pbar + w*(rx*rx');
    Pyy = Pyy + w*(ry*ry');
    Pxy = Pxy + w*(rx*ry');
end

%Measurement update
nu = y - ybar; %Innovation
S = Pyy + R; %Innovation covariance
K = Pxy/S; %Kalman Filter gain

xhat = [xbar(1:3)+K(1:3,:)*nu; qmultiply(gtoq(K(4:6,:)*nu),xbar(4:7)); xbar(8:13)+K(7:12,:)*nu];
Phat = Pbar - K*S*K';

end

% ---------- Helper Functions ---------- %
function q = qmultiply(q1, q2)
    %Quaternion multiplication
    s1 = q1(1);
    v1 = q1(2:4);
    s2 = q2(1);
    v2 = q2(2:4);
    q = [s1*s2 - v1'*v2; s1*v2 + s2*v1 + cross(v1,v2)];
end

function qc = qconj(q)
    %Quaternion conjugate
    qc = [q(1); -q(2:4)];
end

function q = gtoq(r)
    %Converts a Gibbs vector to a quaternion
    s = 1/sqrt(r'*r+1);
    q = [s; s*r];
end

function r = qtog(q)
    %Converts a quaternion to a Gibbs vector
    r = q(2:4)/q(1);
end

function x1 = rkstep(x0,u0,dt)
    %Single step of the 3rd order Runge-Kutta method
    xdot1 = f(x0,u0);
    xdot2 = f(x0+.5*dt*xdot1,u0);
    xdot3 = f(x0-dt*xdot1+2*dt*xdot2,u0);
    x1 = x0 + (dt/6)*(xdot1 + 4*xdot2 + xdot3);
end
