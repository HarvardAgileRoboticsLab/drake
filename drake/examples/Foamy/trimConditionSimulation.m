% Simple script to run TILQR for trim conditions.

function trimConditionSimulation()
    % clear
    clear all; close all; clc;

    % make plane model
    plane = FoamyPlant();

    % find trim condition trajectory
    [xtraj, utraj] = plane.runDircol(0);

    % specify matrices / parameters
    Q = eye(12);
    R = eye(4);
    Qn = Q;
    t0 = xtraj.tspan(1);
    tf = xtraj.tspan(2);
    ts = 0.0001;
    
    % calculate K matrix using Zac's discrete-time code
    %[x0, u0, K_discrete, t, P] = foamy_tvlqr(plane,Q,R,Qn,xtraj,utraj,ts);

    % calculate continuous-time A and B matrices over points along the trajectory
    all_A = [];
    all_B = [];
    for t = linspace(t0, tf, 10)
        [A, B] = calculateAB(plane, t, xtraj.eval(t), utraj.eval(t));
        all_A = [all_A, A(:)];
        all_B = [all_B, B(:)];
    end
    
    % plot A and B values through time -- as a reasonability check, they stay constant
    figure
    subplot(1, 2, 1)
    plot(linspace(t0, tf, 10), all_A')
    title('A values over time')
    subplot(1, 2, 2)
    plot(linspace(t0, tf, 10), all_B')
    title('B values over time')
    
    % run single LQR
    [K, S] = lqr(A, B, Q, R);
    
    %fprintf('Comparison of K (discrete-time TVLQR) and K (continuous-time TILQR):');
    %disp(K_discrete(:, :, 1))
    %disp(K)
    
    % create feedback system with this K matrix
    control = FoamyTILQRController(plane, K, xtraj, utraj, t0, tf);
    tilqr_sys = feedback(plane, control);
    
    % give initial conditions
    x00 = xtraj.eval(0) + [1; 0; 0; zeros([10, 1])];
    
    % simulate!
    xtraj_tilqr = simulate(tilqr_sys, [t0 tf], x00);
    
    % visualize
    t = t0:ts:tf;
    plotting(t, xtraj.eval(t), xtraj_tilqr.eval(t));
    
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
    set(gca, 'xlim', [-6 6], 'ylim', [-6 6], 'zlim', [-6 6]);
    xlabel('x'), ylabel('y'), zlabel('z');
    legend('nominal', 'TILQR', 'Location', 'northeast');
end



function [A, B] = calculateAB(plane, t, x0, u0)
    options.grad_method = 'user_then_numerical';
    [f, df] = geval(@plane.dynamics, t, x0, u0, options);

    sx = numel(x0);

    Aq = df(:, 1 + (1:sx));
    Bq = df(:, sx + 1 + (1:numel(u0)));

    qk = x0(4:7);
    sk = qk(1);
    vk = qk(2:4);

    Gk = [-vk'; sk*eye(3) + hat(vk)];

    % discrete matrices?
    A = blkdiag(eye(3),Gk',eye(6))*Aq*blkdiag(eye(3),Gk,eye(6));
    B = blkdiag(eye(3),Gk',eye(6))*Bq;
end


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

function xhat = hat(x)

xhat = [  0   -x(3)  x(2)
         x(3)   0   -x(1)
        -x(2)  x(1)  0];
end
