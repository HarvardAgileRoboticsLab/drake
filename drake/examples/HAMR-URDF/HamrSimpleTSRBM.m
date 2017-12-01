classdef HamrSimpleTSRBM < TimeSteppingRigidBodyManipulator
    
    properties (SetAccess = protected, GetAccess = public)
        q0
        grav = [0; 0; -9.81e-3];
    end
    
    methods
        
        function obj= HamrSimpleTSRBM(urdf,options)
            
            typecheck(urdf,'char');
            
            obj = obj@TimeSteppingRigidBodyManipulator(urdf,options.dt,options);
            obj.q0 = zeros(obj.getNumPositions(), 1);
            
            %set gravity
            obj.manip = obj.manip.setGravity(obj.grav);
            obj.manip = compile(obj.manip);
            
            % Set initial state
            q0 = zeros(obj.getNumPositions, 1);
            if options.floating
                q0(3) = 12.62;
                q0(8) = 0.520967;
                q0(10) = 0.520967;
                q0(12) = -0.520967;
                q0(14) = -0.520967;
            else
                q0(2) = 0.520967;
                q0(4) = 0.520967;
                q0(6) = -0.520967;
                q0(8) = -0.520967;
            end
            obj.q0 = q0;            
        end        
        
        function obj = compile(obj)
            obj = compile@TimeSteppingRigidBodyManipulator(obj);
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

