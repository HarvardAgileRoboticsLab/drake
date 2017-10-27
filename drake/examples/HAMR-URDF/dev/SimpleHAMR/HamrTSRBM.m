classdef HamrTSRBM < TimeSteppingRigidBodyManipulator
    
    properties (SetAccess = protected, GetAccess = public)
        q0
        grav = [0; 0; -9.81e-3];
    end
    
    methods
        
        function obj=HamrTSRBM(urdf,options)
            
            typecheck(urdf,'char');
            
%             if nargin < 2
%                 options = struct();
%             end
%             if ~isfield(options,'dt')
%                 options.dt = 1;
%             end
%             if ~isfield(options,'floating')
%                 options.floating = true;
%             end
%             if ~isfield(options,'terrain')
%                 options.terrain = RigidBodyFlatTerrain;
%             end


            obj = obj@TimeSteppingRigidBodyManipulator(urdf,options.dt,options);
            obj.q0 = zeros(obj.getNumPositions(), 1);
            
                                    
            %set gravity
            obj.manip = obj.manip.setGravity(obj.grav); 
            obj.manip = compile(obj.manip); 
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

            
%             %Add Ouputs
%             joint_names = obj.getJointNames();
%             actuated_dof = obj.getActuatedJoints();
%             obj = obj.setNumOutputs(obj.getNumStates()+ 2*numel(actuated_dof));
%             state_frame = obj.getStateFrame();
%             
%             if obj.getManipulator().getBody(2).floating
%                 act_jt_names = [joint_names(actuated_dof - 4); ...
%                     joint_names(actuated_dof - 4)];
%                 %                 act_jt_names =
%             else
%                 act_jt_names = [joint_names(actuated_dof + 1); ...
%                     joint_names(actuated_dof +1 )];
%             end
%             
%             output_frame = MultiCoordinateFrame( ...
%                 {state_frame.getFrameByNum(1), ...
%                 state_frame.getFrameByNum(2), ...
%                 CoordinateFrame('ActuatorDeflectionandRate', 2*obj.getNumActuatedDOF(), ...
%                 {}, act_jt_names)},[ones(obj.getNumPositions(),1); ...
%                 2*ones(obj.getNumVelocities(),1); 3*ones(2*numel(actuated_dof),1)]);
%             
%             obj = obj.setOutputFrame(output_frame);
%             
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

%       function varargout = manipulatorDynamics(obj,varargin)
%           q = varargin{1}; v = varargin{2};
%           varargout = cell(1,nargout);
%           [varargout{:}]  = obj.hamr_manip.manipulatorDynamics(q, v);
%           C = varargout{2};
%           varagout{2} = C + obj.loop_forces(q, v);
%       end

%      function Floop = loop_forces(obj, q,v)
%             manip = obj.getManipulator();
%             nBodies = length(manip.body);
%             Floop = zeros(numel(q), 1);
%             net_wrenches = cell(nBodies, 1);
%             for i = 1:nBodies
%                 net_wrenches{i} = zeros(6,1);
%             end
%             %                         Fextg1 = zeros(1, numel(q));
%             %             Fextg2 = zeros(1, numel(q));
%             options.compute_JdotV = true;
%             options.use_mex = false;
%             kinsol = obj.doKinematics(q, v, options);               % kinematics
%             f_ext = zeros(6,getNumBodies(manip));
%             for i = 1:obj.NLOOP
%                 fexti = zeros(6, getNumBodies(manip));
%                 kinopt.base_or_frame_id = obj.LOOP_LINKS(i,1);      % first link in chain
%                 kinopt.rotation_type = 1;                           % we want euler angles
%                 [x1to2, J1to2] = obj.forwardKin(kinsol, obj.LOOP_LINKS(i,2), zeros(3,1), kinopt);
% %                 xd1to2 = J1to2*v;
%
% %                 w1to2 = quatdot2angularvel(x1to2(4:7), xd1to2(4:7));
% %                 aa1to2 = quat2axis(x1to2(4:7));
%
% %                 qinv1to2 = quatConjugate(q1to2);
%
%                 torque = -obj.KD(i)*x1to2(obj.LOOP_AXES(i));  ...
%                     %- obj.BD(i)*xd1to2(obj.LOOP_AXES(i));
%
% %                 torque = -obj.KD(i)*aa1to2(4)*aa1to2(obj.LOOP_AXES(i)) ...
% %                     - obj.BD(i)*w1to2(obj.LOOP_AXES(i));
%
%                 wrench_on_child_in_child_joint_frame = [zeros(2,1);torque;zeros(3,1)];
%
%                 % transform from body frame to joint frame
%                 AdT_body_to_joint = transformAdjoint(manip.body(obj.LOOP_LINKS(i,2)).T_body_to_joint);
%                 fexti(:,obj.LOOP_LINKS(i,2)) = f_ext(:, obj.LOOP_LINKS(i,2)) + ...
%                     AdT_body_to_joint' * wrench_on_child_in_child_joint_frame;
%
%                 if obj.LOOP_LINKS(i,1) ~= 0 % don't apply force to world body
%                     T_parent_body_to_child_joint = homogTransInv(manip.body(obj.LOOP_LINKS(i,2)).Ttree);
%                     AdT_parent_body_to_child_joint = transformAdjoint(T_parent_body_to_child_joint);
%                     fexti(:,obj.LOOP_LINKS(i,1)) = - AdT_parent_body_to_child_joint' * wrench_on_child_in_child_joint_frame;
%                 end
%                 f_ext = f_ext + fexti;
%             end
%             %             f_extN
%             for i = 2:nBodies
%                 external_wrench = f_ext(:, i);
%                 % external wrenches are expressed in body frame. Transform from body to world:
%                 AdT_world_to_body = transformAdjoint(homogTransInv(kinsol.T{i}));
%                 external_wrench = AdT_world_to_body' * external_wrench;
%                 net_wrenches{i} = -external_wrench;
%             end
%
%             for i = nBodies : -1 : 2
%                 body =manip.body(i);
%                 joint_wrench = net_wrenches{i};
%                 Ji = kinsol.J{i};
%                 tau = Ji'*joint_wrench;
%                 Floop(body.velocity_num) = tau;
%                 net_wrenches{body.parent} = net_wrenches{body.parent} + joint_wrench;
%
%             end
% %             Floop
%             %             nWN = reshape(cell2mat(net_wrenches), 6, 9)
%         end