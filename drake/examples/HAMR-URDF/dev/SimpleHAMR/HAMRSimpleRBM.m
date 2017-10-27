classdef HAMRSimpleRBM < RigidBodyManipulator
    
    properties (SetAccess = protected, GetAccess = public)
        q0
        grav = [0; 0; -9.81e-3];
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
                q0(7) = 0.5186;
                q0(9) = 0.5186;
                q0(11) = -0.5186;
                q0(13) = -0.5186;
            else
                q0(1) = 0.5186;
                q0(3) = 0.5186;
                q0(5) = -0.5186;
                q0(7) = -0.5186;
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
