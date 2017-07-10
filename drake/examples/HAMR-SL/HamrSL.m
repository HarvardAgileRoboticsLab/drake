classdef HamrSL < TimeSteppingRigidBodyManipulator
    
    properties (SetAccess = protected, GetAccess = public)
        x0
        grav = [0; 0; -9.81e-3];
    end
    
    methods
        
        function obj=HamrSL(urdf,options)
            
            typecheck(urdf,'char');
            
            if nargin < 2
                options = struct();
            end
            if ~isfield(options,'dt')
                options.dt = 0.001;
            end
            if ~isfield(options,'floating')
                options.floating = true;
            end
            if ~isfield(options,'terrain')
                options.terrain = RigidBodyFlatTerrain;
            end
            
            obj = obj@TimeSteppingRigidBodyManipulator(urdf,options.dt,options);
            obj.x0 = zeros(obj.getNumDiscStates(), 1);
        end
        
        
        function [y, Jy] = output(obj, t, x, u)
            
            nx = numel(x);
            actuated_dof = obj.getActuatedJoints();
            y = [x; x(actuated_dof); x(actuated_dof + nx/2)];
            
            if nargout > 1
                ny = numel(y); nu = numel(u);
                Jy1 = [zeros(nx, 1), eye(nx), zeros(nx, nu)];       % derivative of state
                Jy2 = zeros(ny - nx, 1 + nx + nu);                  % derivative of actuated dof q and qd
                r_ind = (1:ny-nx)';
                c_ind = 1+[actuated_dof; actuated_dof + nx/2];
                Jy2(sub2ind([ny-nx, 1 + nx + nu], r_ind, c_ind)) = 1;
                Jy = [Jy1; Jy2];
            end
        end
        
        function obj = compile(obj)
            obj = compile@TimeSteppingRigidBodyManipulator(obj);
            
            %set gravity
            obj.manip = obj.manip.setGravity(obj.grav);
            
            w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
            obj.manip = compile(obj.manip); 
            warning(w);
            
            %Add Ouputs
            joint_names = obj.getJointNames();
            actuated_dof = obj.getActuatedJoints();
            obj = obj.setNumOutputs(obj.getNumStates()+ 2*numel(actuated_dof));
            state_frame = obj.getStateFrame();
            
            if obj.getManipulator().getBody(2).floating
                act_jt_names = [joint_names(actuated_dof - 4); ...
                    joint_names(actuated_dof - 4)];
                %                 act_jt_names =
            else
                act_jt_names = [joint_names(actuated_dof + 1); ...
                    joint_names(actuated_dof +1 )];
            end
            
            output_frame = MultiCoordinateFrame( ...
                {state_frame.getFrameByNum(1), ...
                state_frame.getFrameByNum(2), ...
                CoordinateFrame('ActuatorDeflectionandRate', 2*obj.getNumActuatedDOF(), ...
                {}, act_jt_names)},[ones(obj.getNumPositions(),1); ...
                2*ones(obj.getNumVelocities(),1); 3*ones(2*numel(actuated_dof),1)]);
            
            obj = obj.setOutputFrame(output_frame);
            
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
            x0 = obj.x0;
        end
        
%         % Added by Neel to change gravity
%         function obj = setGravity(obj, grav)
%             obj.manip = setGravity(obj.manip, grav);
%             obj = compile(obj);
%         end
%         
    end
    
    
    
end

