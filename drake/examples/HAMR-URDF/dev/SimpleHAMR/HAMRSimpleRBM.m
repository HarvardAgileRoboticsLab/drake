classdef HAMRSimpleRBM < RigidBodyManipulator
    
    properties (SetAccess = protected, GetAccess = public)
        q0
        grav = [0; 0; -9.81e-3];
        foot_pos = [0 0 -14.988382167532292];
    end
    
    methods
        
        function obj=HAMRSimpleRBM(urdf,options)
            
            typecheck(urdf,'char');        
            obj = obj@RigidBodyManipulator(urdf, options);
            
            %set gravity
            obj = obj.setGravity(obj.grav); 
            obj = compile(obj); 
            
            % Set initial state
            q0 = zeros(obj.getNumPositions, 1); 
            if options.floating
                q0(3) = 12.5417; 
                q0(8) = 0.520966701443923;
                q0(10) = 0.520966701443923;
                q0(12) = -0.520966701443923;
                q0(14) = -0.520966701443923;
            else
                q0(2) = 0.520966701443923;
                q0(4) = 0.520966701443923;
                q0(6) = -0.520966701443923;
                q0(8) = -0.520966701443923;
            end
            obj.q0 = q0; 
        end   
       
        function nActuatedDOF = getNumActuatedDOF(obj)
            nActuatedDOF = numel(obj.getActuatedJoints());
        end
        
        function x0 = getInitialState(obj)
            q0 = obj.q0;            
            x0 = [q0; 0*q0]; 
        end
        
    end
end
