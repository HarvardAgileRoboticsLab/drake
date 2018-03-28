classdef PendulumTVLQRController < DrakeSystem
    properties
        plant
        xtraj
        utraj
        Q
        R
        Qf
        t0
        tf
        K_array
        S_array
        nx
        nu
        start_time
    end
    methods
        function obj = PendulumTVLQRController(plant, xtraj, utraj, Q, R, Qf, t0, tf)
            
          % call the parent class constructor:
            obj = obj@DrakeSystem(...
             0, ... % number of continuous states
             0, ... % number of discrete states
             26, ... % number of inputs
             4, ... % number of outputs
             true, ... % because the output does not depend on u
             false);  % because the dynamics and output do not depend on t
         
            obj.plant = plant;
            obj.xtraj = xtraj;
            obj.utraj = utraj;
            obj.Q = Q;
            obj.R = R;
            obj.Qf = Qf;
            obj.t0 = t0;
            obj.tf = tf;
            obj.nx = size(xtraj, 1) - 2; %24
            obj.nu = size(utraj, 1); %4
            obj.start_time = t0;
            
            % calculate K array
            [obj.K_array, obj.S_array] = obj.calculate_K_tvlqr(plant, xtraj, utraj, Qf, Q, R);
            
            obj = obj.setInputFrame(plant.getStateFrame);
            obj = obj.setOutputFrame(plant.getInputFrame);
        end
        
        function y = output(obj, tt, ~, x)
            % shift by start_time (hybrid case, but generalizes)
            t = tt - obj.start_time + obj.t0;
            
            if (t > obj.tf)
                warning('Past final time for this controller! before: t=%.5f, after: t=%.5f\n', tt, t);
                t = obj.tf;
            elseif (t < obj.t0)
                warning('Before initial time for this controller! before: t=%.5f, after: t=%.5f\n', tt, t);
                t = obj.t0;
            end
            
            % compute current K, x0, u0
            K = reshape(interp1(obj.K_array.t, obj.K_array.K', t), [obj.nu, obj.nx]);
            x0 = obj.xtraj.eval(t);
            %x0 = convertToFoamyState(x0);
            u0 = obj.utraj.eval(t);
            
            % do quaternion transformation
            q = x(4:7);
            q0 = x0(4:7);
            v1 = [zeros(3,1), eye(3)]*qmultiply(qconj(q0),q);
            
            q = x(11:14);
            q0 = x0(11:14);
            v2 = [zeros(3,1), eye(3)]*qmultiply(qconj(q0),q);

            u = u0 - K*[x(1:3)-x0(1:3); v1; x(8:10)-x0(8:10); v2; x(15:26) - x0(15:26)];
            
            y = applyControlLimits(u);
        end
        
        
        function [K_array_tvlqr, S_array_tvlqr] = calculate_K_tvlqr(obj, plant, xtraj, utraj, Qf, Q, R)

            % dynamics and "initial" conditions for S2
            S0 = reshape(Qf, [obj.nx, obj.nx]);
            
            f = @(t, S) obj.S_dynamics(t, S, plant, Q, R, xtraj, utraj);
            opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-4, 'MaxStep', 0.01);
            S_sol = ode45(f, [obj.t0 obj.tf], S0, opts);

            % flip S back to correct order
            S_sol.y = fliplr(S_sol.y);

            % calculate K at each point
            K_array_tvlqr.t = S_sol.x;
            S_array_tvlqr.t = S_sol.x;
            
            K_array_tvlqr.K = zeros(4 * obj.nx, numel(S_sol.x));
            S_array_tvlqr.S = zeros(obj.nx * obj.nx, numel(S_sol.x));
            
            for i = 1:numel(S_sol.x)
                t = S_sol.x(i);
                x0 = xtraj.eval(t); u0 = utraj.eval(t);
                
                [A, B] = calculate_AB(obj, plant, t, x0, u0);
                
                S = reshape(S_sol.y(:, i), [obj.nx, obj.nx]);
                
                S_array_tvlqr.S(:, i) = reshape(S, [obj.nx * obj.nx, 1]);
                K_array_tvlqr.K(:, i) = reshape(R * B' * S, [obj.nx * 4, 1]);
            end
            
            % plot both S and K
%             figure
%             subplot(1, 2, 1), plot(K_array_tvlqr.t, K_array_tvlqr.K);
%             title('TVLQR: Values in K');
%             subplot(1, 2, 2), plot(S_array_tvlqr.t, S_array_tvlqr.S);
%             title('TVLQR: Values in S');
        end
        
        function dS = S_dynamics(obj, t, S, plant, Q, R, xtraj, utraj)
            x0 = xtraj.eval(t);
            u0 = utraj.eval(t);

            [A, B] = calculate_AB(obj, plant, t, x0, u0);

            S = reshape(S, [obj.nx, obj.nx]);
            dS = Q - S * B * inv(R) * B' * S + S * A + A' * S;
            dS = reshape(dS, [obj.nx * obj.nx, 1]);

            %if (min(eig(S))<0) warning('S is not positive definite, t=%.3f', t); end
        end
        
        function [A, B] = calculate_AB(obj, plant, t, x0, u0)
            options.grad_method = 'user_then_numerical';
            [f, df] = geval(@plant.dynamics, t, x0, u0, options);
            
            sx = size(obj.xtraj, 1);
            
            Aq = df(:, 1 + (1:sx));
            Bq = df(:, sx + 1 + (1:obj.nu));
            
            % first Gk
            qk = x0(4:7);
            sk = qk(1);
            vk = qk(2:4);

            Gk = [-vk'; sk*eye(3) + hat(vk)];
            
            % construct second Gk (for pendulum)
            qp = x0(11:14);
            sp = qp(1);
            vp = qp(2:4);

            Gp = [-vp'; sp*eye(3) + hat(vp)];

            % discrete matrices?
            A = blkdiag(eye(3), Gk', eye(3), Gp', eye(12))*Aq*blkdiag(eye(3), Gk, eye(3), Gp, eye(12));
            B = blkdiag(eye(3), Gk', eye(3), Gp', eye(12))*Bq;
        end
    end
end

function xf = convertToFoamyState(xp)
    xf = [xp(1:7); xp(15:20)];
end


    

% ---------- Helper Functions ---------- %
function y = applyControlLimits(u)
    y(1) = threshold(u(1), 0.00001, 1);
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

