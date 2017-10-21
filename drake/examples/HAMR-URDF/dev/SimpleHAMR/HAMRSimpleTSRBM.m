classdef HAMRSimpleTSRBM < TimeSteppingRigidBodyManipulator
    
    properties (SetAccess = protected, GetAccess = public)
        q0
        grav = [0; 0; -9.81e-3];
    end
    
    methods
        
        function obj=HAMRSimpleTSRBM(urdf,options)
            
            typecheck(urdf,'char');
            
            if nargin < 2
                options = struct();
            end
            if ~isfield(options,'dt')
                options.dt = 1;
            end
            if ~isfield(options,'floating')
                options.floating = true;
            end
            if ~isfield(options,'terrain')
                options.terrain = RigidBodyFlatTerrain;
            end
            
            obj = obj@TimeSteppingRigidBodyManipulator(urdf, options.dt, options);
            
            %set gravity
            obj.manip = obj.manip.setGravity(obj.grav); 
            obj.manip = compile(obj.manip); 
            
            % Set initial state
            q0 = zeros(obj.getNumPositions, 1); 
            if options.floating
                q0(3) = 12.2; 
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
        
        function obj = setInitialState(obj,x0)
            typecheck(x0,'double');
            sizecheck(x0,obj.getNumStates());
            obj.x0 = x0;
        end
        
        function x0 = getInitialState(obj)
            q0 = obj.q0;            
            x0 = [q0; 0*q0]; 
        end
        
    end
end
