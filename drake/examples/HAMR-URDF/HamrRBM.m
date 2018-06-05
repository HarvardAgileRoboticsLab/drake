classdef HamrRBM < RigidBodyManipulator
    
    properties (SetAccess = protected, GetAccess = public)
        x0
        GRAV = [0; 0; -9.81e-3];
        VALID_LOOPS;                                        % index of valid loop constraints
        NL = [];                                            % number of loop constraints
        LEG_NAME = {'FLL4', 'RLL4', 'FRL4', 'RRL4'};        % leg frames
        LEG_ID = [];                                        % leg frame ids
        FOOT_POS = [0, 7.58, -11.350;                       % Position of foot in LEG_NAME frame
            0, 7.58, -11.350;
            0, -7.58, -11.350;
            0, -7.58, -11.350];
        MARKER_POS =[0.06, 8.565, -6.322;                   % Position of vicon marker in LEG_NAME frame
            -0.06, 8.565, -6.322;
            0.06, -8.565, -6.322;
            -0.06, -8.565, -6.322];
        ULIM = 0.3;
        STIFFNESS_MULT = ones(4,1);                         % default (swing act, lift act, swing flex, lift flex)
    end
    
    methods
        
        function obj=HamrRBM(urdf,options)
            
            typecheck(urdf,'char');
            
            obj = obj@RigidBodyManipulator(urdf,options);
            
            if options.floating
                obj.x0 = zeros(2*obj.getNumPositions(), 1);
                obj.x0(3) = 12.69;
            else
                obj.x0 = zeros(2*obj.getNumPositions(), 1);
            end
            
            % stiffness multipliers
            if isfield(options, 'stiffness_mult')
                disp('Modifying stiffness')
                obj.STIFFNESS_MULT = options.stiffness_mult;
                force_objects = obj.force;
                fobj_ind = {25:2:31, 26:2:32, [1:8, 17, 19, 21, 23], [9:16, 18, 20, 22, 24]};   %(swing act, lift act, swing flex, lift flex)

                for i = 1:numel(obj.STIFFNESS_MULT)
                    for j = 1:numel(fobj_ind{i})
                        force_objects{fobj_ind{i}(j)}.k = obj.STIFFNESS_MULT(i)*force_objects{fobj_ind{i}(j)}.k;
                    end
                end
                obj.force = force_objects;
            else
                disp('using default stiffness (specified in URDF)')
            end
            
            % leg frame ids
            for i = 1:numel(obj.LEG_NAME)
                obj.LEG_ID(i) = obj.findLinkId(obj.LEG_NAME{i});
            end
            
            %set gravity
            obj = obj.setGravity(obj.GRAV);
            obj = compile(obj);
            
            % loop const
            valid_loops = [1;2;8;9;13;15];
            obj.VALID_LOOPS = [valid_loops; 18+valid_loops; 36+valid_loops; 54+valid_loops];
            NL = length(obj.VALID_LOOPS); obj.NL = NL;
            
            % set input limits
            umax = obj.ULIM*ones(obj.getNumInputs(), 1);
            umin = -obj.ULIM*ones(obj.getNumInputs(), 1);
            obj = obj.setInputLimits(umin, umax);
            
        end
        
        function obj = compile(obj)
            obj = compile@RigidBodyManipulator(obj);
            
            %Add Ouputs.
            joint_names = {obj.body.jointname}';
            actuated_dof = obj.getActuatedJoints();
            obj = obj.setNumOutputs(obj.getNumStates()+ numel(actuated_dof));
            state_frame = obj.getStateFrame();
            
            if obj.getBody(2).floating
                act_jt_names = joint_names(actuated_dof - 4);
            else
                act_jt_names = joint_names(actuated_dof + 1);
            end
            
            output_frame = MultiCoordinateFrame( ...
                {state_frame.getFrameByNum(1), ...
                state_frame.getFrameByNum(2), ...
                CoordinateFrame('ActuatorDeflection', obj.getNumActuatedDOF(), ...
                {}, act_jt_names)},[ones(obj.getNumPositions(),1); ...
                2*ones(obj.getNumVelocities(),1); 3*ones(numel(actuated_dof),1)]);
            
            obj = obj.setOutputFrame(output_frame);
            
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
        
        function [yfoot, Jfoot, dJfoot] = getFootPosition(obj, q, v, opt)
            
            if strcmpi(opt.loc, 'foot')                     % foot position
                fp = obj.FOOT_POS;
            elseif strcmpi(opt.loc, 'marker')               % marker position
                fp = obj.MARKER_POS;        
            end
            
            fkopt.base_or_frame_id = obj.findLinkId(opt.base_frame);
            
            yfoot = 0*fp;
            Jfoot = cell(size(fp,1), 1);
            dJfoot = cell(size(fp,1), 1);
            
            if nargout > 2
                kinsol = obj.doKinematics(q, v, 'compute_gradients', true);                
                for i = 1:size(fp,1)
                    [yfoot(i,:), Jfoot{i}, dJfoot{i}] = obj.forwardKin(kinsol, obj.LEG_ID(i), fp(i,:)', fkopt);                    
                end
                
            else
                kinsol = obj.doKinematics(q, v);                
                for i = 1:size(fp,1)
                    [yfoot(i,:), Jfoot{i}] = obj.forwardKin(kinsol, obj.LEG_ID(i), fp(i,:)', fkopt);                    
                end
            end        
            yfoot = yfoot';
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
        
    end
end



%% Old Code

%         function [f, df] = dynamics(obj, t, x, u)
%
%             xin = [t; x; u];
%             [f,df] = dynamics_fun(obj,xin);
% %             fprintf('Time: %f \r', t);
%
% %             df_fd = zeros(size(df));
% %             step = max(1e-8, sqrt(eps(max(abs(xin)))));
% %             dxin = step*eye(length(xin));
% %             for k = 1:length(xin)
% %                 df_fd(:,k) = (dynamics_fun(obj, xin+dxin(:,k)) - dynamics_fun(obj, xin-dxin(:,k)))/(2*step);
% %             end
% %
% %             disp('Dynamics Derivative Error:');
% %             disp(max(abs(df_fd(:)-df(:))));
%
%         end
%
%         function[xdot, dxdot] = dynamics_fun(obj,xin)
%
%             % load in parameters
%             nq = obj.getNumPositions();
%             nv = obj.getNumVelocities();
%             nu = obj.nu;
%             nc = obj.nc;
%             nl = obj.nl;
%             nd = obj.nd;
%
%             t = xin(1);
%             q = xin(1+(1:nq));
%             v = xin(1+nq+(1:nv));
%             u = xin(1+nq+nv+(1:nu));
%             l = xin(1+nq+nv+nu+(1:nl));
%             c = xin(1+nq+nv+nu+nl+(1:nc));
%             b = xin(1+nq+nv+nu+nl+nc+(1:nc*nd));
%
%             % Contact basis
%             kinopts = struct();
%             kinopts.compute_gradients = true;
%             kin = obj.doKinematics(q, v, kinopts);
%             [~,~,~,~,~,~,~,~,n,D,dn,dD] = obj.contactConstraints(kin);
%
%             if isempty(n)
%                 n = zeros(0,nq);
%                 dn = zeros(0,nq);
%             end
%
%             D = reshape(cell2mat(D')',nq,nc*nd)';
%             dD = reshape(cell2mat(dD)',nq,nc*nd*nq)';
%
%             % loop constraints
%             [~, K, dK] = obj.positionConstraints(q);
%             K = K(obj.valid_loops, :);
%             dK = reshape(dK(obj.valid_loops, :)', nq, nl*nq)'; %
%
%             % Manipulator Dynamics
%             [H,C,B,dH,dC,dB] = manipulatorDynamics(obj, q, v);
%
%             % form rhs and derivs
%             rhs = B*u + K'*l+ n'*c + D'*b - C;
%
%             drhs_dq = kron(u', eye(nv))*dB(:,1:nq) + kron(c',eye(nq))*comm(nc,nq)*dn + ...
%                 kron(l', eye(nv))*dK + kron(b', eye(nv))*comm(nd*nc,nv)*dD - dC(:,1:nq);
%             drhs_dv = kron(u', eye(nv))*dB(:,nq+(1:nv)) - dC(:, nq+(1:nv));
%
%             % form Hinv and derivs
%             Hinv = inv(H);
%             dH = reshape(dH,nq*nq,nq+nv);
%             dHinv_dq = -kron(Hinv', Hinv)*dH(:,1:nq);
%
%             % Solve for accelerations
%             qddot = Hinv*rhs;
%
%             % xdot and deriv
%             xdot = [v; qddot];
%             dxdot = [zeros(nv, 1), zeros(nv, nq), eye(nv), zeros(nv, nu+nl+nc+nc*nd); ... %dv/dt, dv/dq, dv/dv, dv/du, dv/dl, dv/dc, dv/db
%                 zeros(nv, 1), kron(rhs', eye(nv))*dHinv_dq + Hinv*drhs_dq, Hinv*drhs_dv, Hinv*B, Hinv*K', Hinv*n', Hinv*D']; %da/dt, da/dq, da/dv, da/du, da/dl, da/dc, da/db
%             %
%         end
%
