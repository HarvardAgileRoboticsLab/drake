classdef HAMRPDTracking < DrakeSystem
    
    
    properties (SetAccess=private)
        robot; % to be controlled
        nq;
        nv;
        qa;
        x_des;
        u0;
        Kp;
        Kd;
%         uu;
    end
    
    
    methods
        function obj = HAMRPDTracking(r, u0, x_des, kp, kd)
            % @param r rigid body manipulator instance
            
            input_frame = r.getOutputFrame();
            output_frame = r.getInputFrame();
            
            obj = obj@DrakeSystem(0,0,r.getNumOutputs,r.getNumInputs,true,true);
            
            obj = setInputFrame(obj,input_frame);
            obj = setOutputFrame(obj,output_frame);
            
            obj.robot = r;
            obj.nq = r.getNumPositions();
            obj.nv = r.getNumVelocities();
            obj.qa = r.getActuatedJoints();
            
            obj.Kp = eye(numel(obj.qa))*kp;
            obj.Kd = eye(numel(obj.qa))*kd;            
          
            obj.u0 = u0;            
            obj.x_des = x_des;
            
%             obj.uu = zeros(numel(obj.qa),101); %u0.eval(u0.getBreaks()); 
        end
        
        
        %         function y = output(obj, t, ~, x)
        %
        %             u0 = obj.u0.eval(t);
        %             x_des = obj.x_des.eval(t);
        %             x_hat = x([obj.qa; obj.nq+obj.qa]);
        %
        %             err = x_des - x_hat;
        %
        %             y = u0 + [obj.Kp, obj.Kd]*err;
        %
        % %             disp(
        %
        %             y(y > 0.3) = 0.3;
        %             y(y < -0.3) = -0.3;
        %
        %         end
        
        
        function y = output(obj, t, ~, x)
            
            dim = 3; % 3D
            u0 = obj.u0.eval(t);
            nq = obj.robot.getNumPositions();
            nv = obj.robot.getNumVelocities();
            nc = obj.robot.getNumContactPairs();
            
            x_des = obj.x_des.eval(t);
            q_des = x_des(1:nq);
            qd_des = x_des(nq+(1:nv));
            qleg_des = zeros(dim,  nc);
            vleg_des = qleg_des;
            
            fkopt.base_or_frame_id = obj.robot.findLinkId('Chassis');
            kinsol_des = obj.robot.doKinematics(q_des, qd_des);
            for i = 1:nc
                [qleg_des(:, i), Jleg_des] = obj.robot.forwardKin(kinsol_des, ...
                    obj.robot.findLinkId(obj.robot.legs{i}), obj.robot.pfFull(i,:)', fkopt);
                vleg_des(:,i) = Jleg_des*qd_des;
            end
            
            q = x(1:nq);
            qd = x(nq+(1:nv));
            qleg = zeros(dim, nc);
            vleg = qleg;
            
            kinsol = obj.robot.doKinematics(q, qd);
            for i = 1:nc
                [qleg(:, i), Jleg] = obj.robot.forwardKin(kinsol, ...
                    obj.robot.findLinkId(obj.robot.legs{i}), obj.robot.pfFull(i,:)', fkopt);
                vleg(:,i) = Jleg*qd;
            end            
            
            q_err = qleg_des([1,3], :) - qleg([1,3], :);
            v_err = vleg_des([1,3], :) - vleg([1,3], :);
            fprintf('Leg Error: %f at %f sec \r', max(abs(q_err(:))), t)
            
            Kp = obj.Kp;
            Kp([1, 6, 7, 8], :) = -Kp([1, 6, 7, 8], :);
%             Kp(1:2:end) = 0;
            
            Kd = obj.Kd;
            Kd([1, 6, 7, 8], :) = -Kd([1, 6, 7, 8], :);
%             Kd(1:2:end) = 0;            
            
            y = Kp*q_err(:) + Kd*v_err(:); 
            
            % Input Limits
            y(y > obj.robot.umax) = obj.robot.umax(1);
            y(y < obj.robot.umin) = obj.robot.umin(1);
       
        end
        
        
    end
    
end