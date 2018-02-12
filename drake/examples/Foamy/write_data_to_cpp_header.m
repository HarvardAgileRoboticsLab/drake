%load('K_u_x_data.mat');

[x0, u0, K, t] = foamy_tvlqr(plane,eye(12),eye(4),eye(12),xtraj,utraj,0.05);

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
fileID = fopen('K_header_file_discrete_tvlqr.cpp', 'w');

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