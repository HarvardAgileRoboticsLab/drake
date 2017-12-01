classdef HamrTSRBM < TimeSteppingRigidBodyManipulator
    
    properties (SetAccess = protected, GetAccess = public)
        q0
        grav = [0; 0; -9.81e-3];
        valid_loops; 
    end
    
    methods
        
        function obj=HamrTSRBM(urdf,options)
            
            typecheck(urdf,'char');

            obj = obj@TimeSteppingRigidBodyManipulator(urdf,options.dt,options);
            obj.q0 = zeros(obj.getNumPositions(), 1);
            
                                    
            %set gravity
            obj.manip = obj.manip.setGravity(obj.grav); 
            obj.manip = compile(obj.manip); 
            
            % set loops
            valid_loops = [1;2;8;9;13;15];
            obj.valid_loops = [valid_loops; 18+valid_loops; 36+valid_loops; 54+valid_loops];  

        end
        
        
%         function [y, Jy] = output(obj, t, x, u)
%             
%             nx = numel(x);
%             actuated_dof = obj.getActuatedJoints();
%             y = [x; x(actuated_dof); x(actuated_dof + nx/2)];
%             
%             if nargout > 1
%                 ny = numel(y); nu = numel(u);
%                 Jy1 = [zeros(nx, 1), eye(nx), zeros(nx, nu)];       % derivative of state
%                 Jy2 = zeros(ny - nx, 1 + nx + nu);                  % derivative of actuated dof q and qd
%                 r_ind = (1:ny-nx)';
%                 c_ind = 1+[actuated_dof; actuated_dof + nx/2];
%                 Jy2(sub2ind([ny-nx, 1 + nx + nu], r_ind, c_ind)) = 1;
%                 Jy = [Jy1; Jy2];
%             end
%         end
        
        function obj = compile(obj)
            obj = compile@TimeSteppingRigidBodyManipulator(obj);

%             
%             %Add Ouputs.
%             joint_names = obj.getJointNames();
%             actuated_dof = obj.getActuatedJoints();
%             obj = obj.setNumOutputs(obj.getNumStates()+ numel(actuated_dof));
%             state_frame = obj.getStateFrame();
%             
%             if obj.getManipulator().getBody(2).floating
%                 act_jt_names = joint_names(actuated_dof - 4);                
%             else
%                 act_jt_names = joint_names(actuated_dof + 1);
%             end
%             
%             output_frame = MultiCoordinateFrame( ...
%                 {state_frame.getFrameByNum(1), ...
%                 state_frame.getFrameByNum(2), ...
%                 CoordinateFrame('ActuatorDeflection', obj.getNumActuatedDOF(), ...
%                 {}, act_jt_names)},[ones(obj.getNumPositions(),1); ...
%                 2*ones(obj.getNumVelocities(),1); 3*ones(numel(actuated_dof),1)]);
%             
%             obj = obj.setOutputFrame(output_frame);
%             
        end       
        
        function hamr_pd = pdtrajtracking(obj, kp, kd)
            
            qa = obj.getActuatedJoints();
            
            Kp = kp*eye(numel(qa)); 
            Kd = kd*eye(numel(qa)); 
            
            hamr_pd = obj.pdcontrol(Kp, Kd, qa);             
            
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
            x0 = [obj.q0; 0*obj.q0]; 
        end

    end    
end

