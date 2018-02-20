classdef FoamyDiscreteTVLQRController < DrakeSystem
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
        function obj = FoamyDiscreteTVLQRController(plant, xtraj, utraj, Q, R, Qf, t0, tf)
            
          % call the parent class constructor:
            obj = obj@DrakeSystem(...
             0, ... % number of continuous states
             0, ... % number of discrete states
             13, ... % number of inputs
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
            obj.nx = size(xtraj, 1)-1;
            obj.nu = size(utraj, 1);
            obj.start_time = t0;
            
            % calculate K array
            ts = (tf - t0)/400;
            [~, ~, KK, Kt, SS] = foamy_tvlqr(plant, Q * ts, R * ts, Qf, xtraj, utraj, ts);
            
            % reshape into consistent form, match dimensions with time
            obj.K_array.K = reshape(KK, [obj.nx * obj.nu, numel(Kt)-1]);
            obj.K_array.K = [obj.K_array.K(:, 1), obj.K_array.K];
            obj.K_array.t = Kt;
            obj.S_array.S = reshape(SS, [obj.nx * obj.nx, numel(Kt)]);
            obj.S_array.t = Kt;
            
            % frames
            obj = obj.setInputFrame(plant.getStateFrame);
            obj = obj.setOutputFrame(plant.getInputFrame);
         
            % plot both S and K
            figure
            subplot(1, 2, 1), plot(obj.K_array.t, obj.K_array.K);
            title('Discrete TVLQR: Values in K');
            subplot(1, 2, 2), plot(obj.S_array.t, obj.S_array.S);
            title('Discrete TVLQR: Values in S');
        end
    
        function y = output(obj, t, ~, x)
            % shift by start_time (hybrid case, but generalizes)
            t = t - obj.start_time + obj.t0;
            
            if (t > obj.tf)
                warning('Past final time for this controller!');
                t = obj.tf;
            elseif (t < obj.t0)
                warning('Before initial time for this controller!');
                t = obj.t0;
            end

            % compute current K, x0, u0
            K = reshape(interp1(obj.K_array.t, obj.K_array.K', t), [obj.nu, obj.nx]);
            x0 = obj.xtraj.eval(t);
            u0 = obj.utraj.eval(t);

            % do quaternion transformation
            q = x(4:7);
            q0 = x0(4:7);
            v = [zeros(3,1), eye(3)]*qmultiply(qconj(q0),q);

            u = u0 - K*[x(1:3)-x0(1:3); v; x(8:13)-x0(8:13)]; 

            y = applyControlLimits(u);
            
            if (t > 1.9942838)
                m = 2;
            end
            
            
        end
    end
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
