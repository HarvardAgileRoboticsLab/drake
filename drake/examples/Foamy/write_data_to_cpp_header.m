%load('K_u_x_data.mat');
tvlqr = 0;

p = FoamyPlant();
[xtraj,utraj] = p.runDircol_trim2();
num_Ks = 100;
ts = (xtraj.tspan(2)-xtraj.tspan(1))/num_Ks;
%ts = .05;
%Q = (1/.05)^2*eye(12);
weights =...
    [0.05,0.05,0.05,... %Position Weights
    10,10,10,...%Orientation Weights
    0.1,0.1,0.1,...%Linear velocity weights
    1,1,1];%Angular Velocity Weights
Q = diag(weights);
Qn = (1/.01)^2*eye(12);
R = (1/.25)^2*eye(4);
%R = 100*eye(4);
if(tvlqr)
    [x0,u0,K,t] = foamy_tvlqr(p,Q,R,Qn,xtraj,utraj,ts);

else
    

    t = xtraj.tspan(1):ts:xtraj.tspan(2);
    
    N = length(t);
    x0 = xtraj.eval(t);
    u0 = utraj.eval(t);
    A = zeros(12,12,N);
    B = zeros(12,4,N);
    K = zeros(4,12,N-1);
    for k = 1:(N-1)

        [~,dxx] = p.dynamics(10,x0(:,k),u0(:,k));
         Aq = dxx(:,2:14);
         Bq = dxx(:,15:18);
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

        [K(:,:,k)] = lqr(A(:,:,k),B(:,:,k),Q,R);
    end
end

%[x0, u0, K, t] = foamy_tvlqr(plane,eye(12),eye(4),eye(12),xtraj,utraj,0.05);

x0 = x0(:, 1:(end-1));
u0 = u0(:, 1:(end-1));

nx = size(x0, 1);
m = size(u0, 1);

nk = size(K, 2);
N = size(x0, 2);

% check dims
if (N ~= size(u0, 2) || N ~= size(K, 3) || m ~= size(K, 1))
    warning('Dimension mismatch?');
end

K = reshape(K, [nk * m, N]);

% open new file
fileID = fopen('K_header_file_discrete_tvlqr_columns.cpp', 'w');

% write t
fprintf(fileID, 'double t[%d] = {', N);
for i = 1:N-1
    fprintf(fileID, '%.5f, ', t(i));
end
fprintf(fileID, '%.5f};\n\n ', t(N));

% write x
fprintf(fileID, 'double x0[%d][%d] = {', nx, N);
for i = 1:nx-1
    fprintf(fileID, '{');
    for j = 1:N-1
        fprintf(fileID, '%.5f, ', x0(i, j));
    end
    fprintf(fileID, '%.5f},\n', x0(i, N));
end
fprintf(fileID, '{');
for j = 1:N-1
    fprintf(fileID, '%.5f, ', x0(nx, j));
end
fprintf(fileID, '%.5f}}; \n\n', x0(nx, N));

% write u
fprintf(fileID, 'double u0[%d][%d] = {', m, N);
for i = 1:m-1
    fprintf(fileID, '{');
    for j = 1:N-1
        fprintf(fileID, '%.5f, ', u0(i, j));
    end
    fprintf(fileID, '%.5f},\n', u0(i, N));
end
fprintf(fileID, '{');
for j = 1:N-1
    fprintf(fileID, '%.5f, ', u0(m, j));
end
fprintf(fileID, '%.5f}}; \n\n', u0(m, N));

% write K array
fprintf(fileID, 'double K[%d][%d] = {', nk * m, N);
for i = 1:(nk * m)-1
    fprintf(fileID, '{');
    for j = 1:N-1
        fprintf(fileID, '%.5f, ', K(i, j));
    end
    fprintf(fileID, '%.5f},\n', K(i, N));
end
fprintf(fileID, '{');
for j = 1:N-1
    fprintf(fileID, '%.5f, ', K(nk * m, j));
end
fprintf(fileID, '%.5f}}; \n\n', K(nk * m, N));

fclose(fileID);

function xhat = hat(x)

xhat = [  0   -x(3)  x(2)
         x(3)   0   -x(1)
        -x(2)  x(1)  0];
end