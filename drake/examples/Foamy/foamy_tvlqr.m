function [x0, u0, K, t, P] = foamy_tvlqr(plant,Q,R,Qn,xtraj,utraj,ts)

t = xtraj.tspan(1):ts:xtraj.tspan(2);
N = length(t);

x0 = xtraj.eval(t);
u0 = utraj.eval(t);

%Calculate discrete-time A and B matrices
A = zeros(12,12,N);
B = zeros(12,4,N);
for k = 1:(N-1)
    [~,Aq,Bq] = rkstep(x0(:,k),u0(:,k),ts);
    
    qk = x0(4:7,k);
    sk = qk(1);
    vk = qk(2:4);
    
    qn = x0(4:7,k+1);
    sn = qn(1);
    vn = qn(2:4);
    
    Gk = [-vk'; sk*eye(3) + hat(vk)];
    Gn = [-vn'; sn*eye(3) + hat(vn)];
    
    A(:,:,k) = blkdiag(eye(3),Gn',eye(6))*Aq*blkdiag(eye(3),Gk,eye(6));
    B(:,:,k) = blkdiag(eye(3),Gn',eye(6))*Bq;
end

%Calculate LQR gain matrices
P = zeros(12,12,N);
K = zeros(4,12,N-1);
P(:,:,N) = Qn;
for k = (N-1):-1:1
    K(:,:,k) = (R + B(:,:,k)'*P(:,:,k+1)*B(:,:,k))\(B(:,:,k)'*P(:,:,k+1)*A(:,:,k));
    P(:,:,k) = Q + K(:,:,k)'*R*K(:,:,k) + (A(:,:,k)-B(:,:,k)*K(:,:,k))'*P(:,:,k+1)*(A(:,:,k)-B(:,:,k)*K(:,:,k));
end

function [xn, A, B] = rkstep(x,u,dt)
    %Single step of the 3rd order Runge-Kutta method with derivatives
    
    [xdot1, dxdot1] = plant.dynamics(0,x,u);
    [xdot2, dxdot2] = plant.dynamics(0,x+.5*dt*xdot1,u);
    [xdot3, dxdot3] = plant.dynamics(0,x-dt*xdot1+2*dt*xdot2,u);
    xn = x + (dt/6)*(xdot1 + 4*xdot2 + xdot3);
    
    A1 = dxdot1(:,2:14);
    B1 = dxdot1(:,15:18);
    A2 = dxdot2(:,2:14);
    B2 = dxdot2(:,15:18);
    A3 = dxdot3(:,2:14);
    B3 = dxdot3(:,15:18);
    
    A = eye(13) + (dt/6)*A1 + (2*dt/3)*A2*(eye(13)+(dt/2)*A1)...
                + (dt/6)*A3*(eye(13) - dt*A1 + 2*dt*A2*(eye(13)+(dt/2)*A1));
    
    B = (dt/6)*B1 + (2*dt/3)*(B2 + A2*(dt/2)*B1)...
                  + (dt/6)*(B3 + A3*(2*dt*(B2 + A2*(dt/2)*B1) - dt*B1));
end

function xhat = hat(x)

xhat = [  0   -x(3)  x(2)
         x(3)   0   -x(1)
        -x(2)  x(1)  0];
end

end


