function write_data_to_cpp_header()
    clear all; close all;

    %% ~~~~~ Generate data ~~~~~

    % flag to represent whether to run TVLQR, or just calculate an array of
    % TILQR K matrices at each time point
    tvlqr = 0;

    % initialize plane model and run trajectories
    p = FoamyPlant();
    [xtraj, utraj] = p.runDircol_trim2();

    % parameters
    num_Ks = 100;
    ts = (xtraj.tspan(2)-xtraj.tspan(1))/num_Ks;

    %ts = .05;
    %Q = (1/.05)^2*eye(12);

    % cost matrices
    weights =...
        [0.05,0.05,0.05,... %Position Weights
        10,10,10,...%Orientation Weights
        0.1,0.1,0.1,...%Linear velocity weights
        1,1,1];%Angular Velocity Weights
    Q = diag(weights);
    %Qn = (1/.01)^2*eye(12);
    %R = (1/.25)^2*eye(4);
    R = 400*eye(4);

    if (tvlqr)
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
            % calculate discrete A and B at each time step
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

            % convert to 12-dim from 13-dim quaternion representation
            A(:,:,k) = blkdiag(eye(3),Gn',eye(6))*Aq*blkdiag(eye(3),Gk,eye(6));
            B(:,:,k) = blkdiag(eye(3),Gn',eye(6))*Bq;

            % stack into array
            [K(:,:,k)] = lqr(A(:,:,k), B(:,:,k), Q, R);
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

    %% ~~~~~ Save data ~~~~~
    
    % save current date and time so that .mat and .cpp match
    current_date_time = datestr(now, 'mm_dd-HH:MM:SS');

    % open new file
    cpp_filename = sprintf('K_header_file_%s.cpp', current_date_time);
    matlab_filename = sprintf('K_header_file_%s.mat', current_date_time);
    
    % write data to header file and matlab file
    write_data = 1;
    
    if (write_data)
        write_to_cpp(cpp_filename, N, nx, nk, m, t, x0, u0, K, Q, R);
        save(matlab_filename);
    end
    
    %% ~~~~~ Simulate data ~~~~~
    
    % next, try to simulate the trajectory forward without Drake, subject
    % to these matrices
    % set initial condition
    x00 = x0(:, 1); x00(1:3) = x00(1:3) + [0.1; 0.01; 0.01];

    % set integration options
    options = odeset('RelTol', 1e-3, 'AbsTol', 1e-3, 'MaxStep', 0.5e-3);

    % simulate forward
    [t_sol, x_sol] = ode45(@(t_curr, x_curr) dynamics(t_curr, x_curr, t, x0, u0, K, p), ...
        [t(1), t(end)], x00, options);
    
    % plot and visualize
    vis = FoamyVisualizer(p);
    xtraj_tilqr = PPTrajectory(spline(linspace(t_sol(1), t_sol(end), numel(t_sol)), x_sol'));
    xtraj_tilqr = xtraj_tilqr.setOutputFrame(vis.getInputFrame);
    vis.playback(xtraj_tilqr);
    
    % plot trajectories
    plotting(t, xtraj.eval(t), xtraj_tilqr.eval(t));
    
    % pause!
    keyboard;
end

function plotting(t, xtraj_eval, xtraj_tilqr_eval)
    % define states for labels
    states = {'x', 'y', 'z', 'q1', 'q2', 'q3', 'q4', 'vx', 'vy', 'vz', 'w1', 'w2', 'w3'};

    % plot all states for the nominal and controlled trajectory
    figure('Position', [100 100 1000 800]);
    for i = 1:12
        subplot(4, 4, i)
        plot(t, xtraj_eval(i,:), 'r'), hold on
        plot(t, xtraj_tilqr_eval(i,:), 'k'),
        title(states{i}, 'Fontsize', 14);
        set(gca, 'Xlim', [t(1), t(end)]);
        ytickformat('%.2f');
    end
    subplot(4, 4, [14 15])
    plot(t, xtraj_eval(i,:), 'r'), hold on
    plot(t, xtraj_tilqr_eval(i,:), 'k'),
    title(states{13}, 'Fontsize', 14);
    set(gca, 'Xlim', [t(1), t(end)]);
    legend('nominal', 'TILQR', 'Location', 'bestoutside');
    ytickformat('%.2f');
    
    % plot 3D trajectory in space
    figure
    plot3(xtraj_eval(1,:), xtraj_eval(2,:), xtraj_eval(3,:), 'r'); hold on;
    plot3(xtraj_tilqr_eval(1,:), xtraj_tilqr_eval(2,:), xtraj_tilqr_eval(3,:), 'k'), hold off
    title('Plane trajectory in space');
    %set(gca, 'xlim', [-6 6], 'ylim', [-6 6], 'zlim', [-6 6]);
    xlabel('x'), ylabel('y'), zlabel('z');
    legend('nominal', 'TILQR', 'Location', 'northeast');
end

function dx = dynamics(t, x, t_array, xd_array, ud_array, K_array, plane)
    xd = interp1(t_array(1:end-1), xd_array', t)';
    ud = interp1(t_array(1:end-1), ud_array', t)';
    K = interp1(t_array(1:end-1), K_array', t);
    
    K = reshape(K, [4, 12]);
    u_control = foamy_controller(x, xd, ud, K);
    u_applied = applyControlLimits(u_control);
    
    dx = plane.dynamics(t, x, u_applied);
end

function y = applyControlLimits(u)
    y(1) = threshold(u(1), 0, 1);
    y(2) = threshold(u(2), -1, 1);
    y(3) = threshold(u(3), -1, 1);
    y(4) = threshold(u(4), -1, 1);
    y = y';
end

function y = threshold(yy, ymin, ymax)
    if (yy < ymin)
        y = ymin;
    elseif (yy > ymax)
        y = ymax;
    else
        y = yy;
    end
end

function write_to_cpp(filename, N, nx, nk, m, t, x0, u0, K, Q, R)
    % open file
    fileID = fopen(filename, 'w');
    
    fprintf(fileID, '// This header file was generated by the "write_data_to_cpp_header.m" function.\n');
    fprintf(fileID, '// --------\n');
    fprintf(fileID, '// K was calculated with the following cost matrices:\n');
    fprintf(fileID, '// Q = diag(%.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f)\n', diag(Q));
    fprintf(fileID, '// R = diag(%.3f, %.3f, %.3f, %.3f)\n', diag(R));
    fprintf(fileID, '\n\n');

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

    % close file
    fclose(fileID);
end


function xhat = hat(x)

xhat = [  0   -x(3)  x(2)
         x(3)   0   -x(1)
        -x(2)  x(1)  0];
end