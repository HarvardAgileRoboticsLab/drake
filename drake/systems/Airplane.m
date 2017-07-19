classdef Airplane < DrakeSystem
    %AIRPLANE Is a base class for implementing fixed-wing aircraft plants
    %   Subclasses must implement an aerodynamics method that computes
    %   forces and torques.
    
    properties
        g = 9.81; %gravitational acceleration (m/s^2)
        rho = 1.225; %Air density (kg/m^3)
    end
    
    methods
        function obj = Airplane()
            obj = obj@DrakeSystem(13,0,4,13,0,1);
            obj = obj.setStateFrame(CoordinateFrame('AirplaneState',13,'x',{'x1','x2','x3','q0','q1','q2','q3','v1','v2','v3','w1','w2','w3'}));
            obj = obj.setInputFrame(CoordinateFrame('AirplaneInput',4,'u',{'thr','ail','elev','rud'}));
            obj = obj.setOutputFrame(obj.getStateFrame); %full state feedback
            %quaternion unit-norm constraint
            obj = obj.addStateConstraint(QuadraticConstraint(.5,.5,blkdiag(zeros(3),eye(4),zeros(6)),zeros(13,1)));
        end
        
        function [xdot,dxdot] = dynamics(obj,t,x,u)
            %States
            %pos = x(1:3); %position in world frame
            q = x(4:7); %rotation from body to world frame, [s; v]
            v = x(8:10); %velocity in world frame
            w = x(11:13); %angular velcity in body frame
            
            %Inputs
            thr = u(1); %Thrust
            ail = u(2); %Aleron
            ele = u(3); %Elevator
            rud = u(4); %Rudder
            
            %Compute aerodynamic forces + torques
            [F,T,dF,dT] = aerodynamics(obj,t,q,v,w,thr,ail,ele,rud);
            
            %Kinematics
            xdot(1:3) = v;
            xdot(4:7) = .5*[-q(2:4)'*w; q(1)*w + cross(q(2:4),w)];
            xdot(8:10) = obj.qrotate(q,F)./obj.mass - [0; 0; obj.g];
            xdot(11:13) = obj.Jinv*(T - cross(w,obj.J*w));
            
            %Derivatives
            dxdot = zeros(13,17);
            dxdot(1:3,8:10) = eye(3);
            dxdot(4:7,4:7) = .5*[0 -w'; w, -obj.hat(w)];
            dxdot(4:7,11:13) = .5*[-q(2:4)'; q(1)*eye(3) + obj.hat(q(2:4))];
            dxdot(8:10,:) = [zeros(3,3), (1/obj.mass)*obj.dqrotate(q,F), zeros(3,10)] + [zeros(3,3), obj.qtodcm(q)*dF];
            dxdot(11:13,:) = [zeros(3,10), obj.Jinv*(obj.hat(obj.J*w) - obj.hat(w)*obj.J)] + [zeros(3,3), obj.Jinv*dT];
        end
        
        function [rrot] = qrotate(obj,q,r)
            %Rotate vector r by quaternion q
            rrot = r + 2*cross(q(2:4),(cross(q(2:4),r) + q(1)*r));
        end
        
        function drdq = dqrotate(obj,q,r)
            drdq = 2*[cross(q(2:4),r); obj.hat(cross(r,q(2:4)))-obj.hat(v)*obj.hat(r)];
        end
        
        function qc = qconj(obj,q)
            %Quaternion conjugate
            qc = [q(1); -q(2:4)];
        end
        
        function R = qtodcm(obj,q)
            %Quaternion to rotation matrix
            R = eye(3) + 2*obj.hat(q(2:4))*(obj.hat(q(2:4)) + eye(3)*q(1));
        end
        
        function xhat = hat(obj,x)
            %The hat map - a 3x3 skew-symmetric matrix equivalent to the
            %cross-product
            xhat = [  0   -x(3)  x(2)
                     x(3)   0   -x(1)
                    -x(2)  x(1)  0];
        end
    end
    
    methods (Abstract=true)
        [F,T,dF,dT] = aerodynamics(obj,t,q,v,w,thr,ail,ele,rud)
    end
    
end

