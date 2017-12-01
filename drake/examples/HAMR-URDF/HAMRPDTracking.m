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
    end
    
    
    methods
        function obj = HAMRPDTracking(r, u0, x_des, kp, kd)
            % @param r rigid body manipulator instance
            
            input_frame = r.getStateFrame();
            output_frame = r.getInputFrame();
            
            obj = obj@DrakeSystem(0,0,r.getNumStates,r.getNumInputs,true,true);
            
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
        end
        
        
        function y = output(obj, t, ~, x)
            
            u0 = obj.u0.eval(t);
            x_des = obj.x_des.eval(t);
            x_hat = x([obj.qa; obj.nq+obj.qa]);
            
            err = x_des - x_hat;
            
            y = u0 + [obj.Kp, obj.Kd]*err;
            
%             disp(
            
            y(y > 0.3) = 0.3;
            y(y < -0.3) = -0.3;             
            
        end
    end
    
end