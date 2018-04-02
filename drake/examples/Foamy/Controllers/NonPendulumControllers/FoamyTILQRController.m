classdef FoamyTILQRController < DrakeSystem
    properties
        plant
        xtraj
        utraj
        K
        t0
        tf
        nx
        nu
    end
    methods
        function obj = FoamyTILQRController(plant, K, xtraj, utraj, t0, tf)
            
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
            obj.K = K;
            obj.t0 = t0;
            obj.tf = tf;
            obj.nx = 12; %size(xtraj, 1)-1;
            obj.nu = 4; %size(utraj, 1);
            
            obj = obj.setInputFrame(plant.getStateFrame);
            obj = obj.setOutputFrame(plant.getInputFrame);
        end
        
        function y = output(obj, t, ~, x)            
            if (t > obj.tf)
                warning('Past final time for this controller! before: t=%.5f, after: t=%.5f\n', tt, t);
                t = obj.tf;
            elseif (t < obj.t0)
                warning('Before initial time for this controller! before: t=%.5f, after: t=%.5f\n', tt, t);
                t = obj.t0;
            end
            
            % compute current x0, u0
            x0 = obj.xtraj.eval(t);
            u0 = obj.utraj.eval(t);
            
            % do quaternion transformation
            q = x(4:7);
            q0 = x0(4:7);
            v = [zeros(3,1), eye(3)]*qmultiply(qconj(q0),q);

            u = u0 - obj.K * [x(1:3)-x0(1:3); v; x(8:13)-x0(8:13)]; 
            
            y = applyControlLimits(u);
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

