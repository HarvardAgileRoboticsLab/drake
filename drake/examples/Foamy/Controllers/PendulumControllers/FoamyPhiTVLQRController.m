classdef FoamyPhiTVLQRController < DrakeSystem
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
        switch_position
        phase
    end
    methods
        function obj = FoamyPhiTVLQRController(plant, xtraj, utraj, Q, R, Qf, t0, tf)
            
          % call the parent class constructor:
            obj = obj@DrakeSystem(...
             0, ... % number of continuous states
             0, ... % number of discrete states
             13, ... % number of inputs
             4, ... % number of outputs
             true, ... % because the output does not depend on u
             false);  % because the dynamics and output do not depend on t
            
            obj = obj.setInputFrame(plant.getStateFrame);
            obj = obj.setOutputFrame(plant.getInputFrame);
         
            obj.plant = plant;
            obj.xtraj = xtraj;
            obj.utraj = utraj;
            obj.Q = Q;
            obj.R = R;
            obj.Qf = Qf;
            obj.t0 = t0;
            obj.tf = tf;
            obj.nx = size(Q, 1);
            obj.nu = size(R, 1);
            
            % make the last point the impact point
            last_point = xtraj.eval(tf);
            obj.switch_position = last_point(1:3);

            % define time array and calculate phi
            t = linspace(t0, tf, 200);
            phi = zeros([1, numel(t)]);
            for i = 1:numel(t)
                phi(i) = obj.calculate_phi(xtraj.eval(t(i)));
            end

            % find dt/dphi and store everything in "phase" structure
            phase.phi = phi;
            phase.t = t;
            phase.dt_dphi = derivatives(t, phi);
            obj.phase = phase;
            
            % plot
            figure
            subplot(1, 2, 1), plot(phase.t, phase.phi, '-*'), xlabel('t'), ylabel('phi')
            subplot(1, 2, 2), plot(phase.phi, phase.dt_dphi, '-*'), xlabel('phi'), ylabel('dt_dphi')
            
            % calculate K array
            [obj.K_array, obj.S_array] = obj.calculate_K_phi(plant, xtraj, utraj, Qf, Q, R);
        end
        
        function y = output(obj, t, ~, x)
            % calculate current phi
            phi = obj.calculate_phi(x);
            
            % get desired x and u
            if (phi < min(obj.phase.phi))
                equivalent_t = obj.t0;
                phi = min(obj.phase.phi);
            else
                equivalent_t = interp1(obj.phase.phi, obj.phase.t, phi);
            end
            xd = obj.xtraj.eval(equivalent_t);
            ud = obj.utraj.eval(equivalent_t);
            current_K = reshape(interp1(obj.K_array.phi, obj.K_array.K', phi), [obj.nu, obj.nx]);
            
            % do quaternion transformation
            q = x(4:7);
            q0 = xd(4:7);
            v = [zeros(3,1), eye(3)] * qmultiply(qconj(q0),q);

            % control
            u = ud - current_K * [x(1:3)-xd(1:3); v; x(8:13)-xd(8:13)];
            
            y = applyControlLimits(u);
            %disp(u')
            
            %if (sum(isnan(y)) > 0)
            %    m = 2;
            %end
        end
        
        function [K_array_phi, S_array] = calculate_K_phi(obj, plant, xtraj, utraj, Qf, Q, R)

            % dynamics and "initial" conditions for S2
            S0 = reshape(Qf, [obj.nx * obj.nx, 1]);
            f = @(phi, S) obj.S_dynamics(phi, S, plant, Q, R, xtraj, utraj);
            opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-4, 'MaxStep', 0.01);
            S_sol = ode45(f, [obj.phase.phi(1) obj.phase.phi(end)], S0, opts);

            % flip S back to correct order
            S_sol.y = fliplr(S_sol.y);

            figure
            plot(S_sol.x, S_sol.y);
            title('Values in S');

            % calculate K at each point
            S_array.phi = S_sol.x;
            S_array.S = S_sol.y;
            
            K_array_phi.phi = S_sol.x;
            K_array_phi.K = zeros(obj.nu * obj.nx, numel(S_sol.x));
            
            for i = 1:numel(S_sol.x)
                % find A and B matrices to compute K
                phi = S_sol.x(i);
                equivalent_t = interp1(obj.phase.phi, obj.phase.t, phi);
                x0 = xtraj.eval(equivalent_t);
                u0 = utraj.eval(equivalent_t);
                
                [A, B] = calculate_AB(plant, equivalent_t, x0, u0);
                
                S = reshape(S_sol.y(:, i), [obj.nx, obj.nx]);

                K_array_phi.K(:, i) = reshape(R * B' * S, [obj.nx * obj.nu, 1]);
                
                K_array_phi.eq_t(:, i) = equivalent_t;
            end
            
            figure
            plot(K_array_phi.phi, K_array_phi.K);
            title('Values in K');
        end
        
        function dS = S_dynamics(obj, phi, S, plant, Q, R, xtraj, utraj)
            equivalent_t = interp1(obj.phase.phi, obj.phase.t, phi);
            dt_dphi = interp1(obj.phase.phi, obj.phase.dt_dphi, phi);
            
            x0 = xtraj.eval(equivalent_t);
            u0 = utraj.eval(equivalent_t);
            
            [A, B] = calculate_AB(plant, equivalent_t, x0, u0);

            S = reshape(S, [obj.nx obj.nx]);

            dS = Q - S * B * inv(R) * B' * S + S * A + A' * S;
            dS = reshape(dS, [obj.nx * obj.nx, 1]);
            dS = dt_dphi * dS; % scaling! this is the phase-based part.

            if (min(eig(S))<0)
                warning('S is not positive definite');
            end
        end
            
        
        % main phase calculatin method (for consistency)
        function phi = calculate_phi(obj, x)
            dist = x(1:3) - obj.switch_position;
            phi = sum(dist.^2);
            phi = -sqrt(phi);
        end
    
    end
end

function [A, B] = calculate_AB(plant, t, x0, u0)
    options.grad_method = 'user_then_numerical';
    [f, df] = geval(@plant.dynamics, t, x0, u0, options);

    sx = numel(x0);
    nu = numel(u0);

    Aq = df(:, 1 + (1:sx));
    Bq = df(:, sx + 1 + (1:nu));

    qk = x0(4:7);
    sk = qk(1);
    vk = qk(2:4);
    
    Gk = [-vk'; sk*eye(3) + hat(vk)];

    % discrete matrices?
    A = blkdiag(eye(3),Gk',eye(6)) * Aq * blkdiag(eye(3),Gk,eye(6));
    B = blkdiag(eye(3),Gk',eye(6)) * Bq;
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
