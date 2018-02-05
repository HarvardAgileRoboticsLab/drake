classdef HamrRBM < RigidBodyManipulator
    
    properties (SetAccess = protected, GetAccess = public)
        q0
        grav = [0; 0; -9.81e-3];    
        pf = [0, 7.58, -11.350;     % Position of foot in local frame (FL; RL; FR; RR)
            0, 7.58, -11.350;
            0, -7.58, -11.35
            0, -7.58, -11.350];        
        foot = {'FLL4', 'RLL4', 'FRL4', 'RRL4'};          % name of foot link
        nu              % number of actual inputs
        loop_const      % loop constraints
        nl              % number of loop constraints
        nc              % number of contact pairs
        nd              % number of basis vectors in polyhedral friction cone
    end
    
    methods
        
        function obj= HamrRBM(options)
            
            urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', 'HAMR_scaledV2.urdf');
            obj = obj@RigidBodyManipulator(urdf,options);
            
            obj.nu = obj.getNumInputs();            
            
            % Change gravity
            obj = obj.setGravity(obj.grav);
            obj = compile(obj);
            
            % contact forces
            [phi,~,d] = obj.contactConstraints(getZeroConfiguration(obj));
            nc = length(phi); obj.nc = nc;
            nd = 2*length(d); obj.nd = nd;
            
            % loop const
            [k, ~] = obj.positionConstraints(getZeroConfiguration(obj));
            nl = length(k); obj.nl = nl;
            obj.loop_const = obj.position_constraints;
            
            % Add contact and loop const as inputs
            obj = obj.setNumInputs(obj.nu + nc + nc*nd + nl);
            
            % Initial state
            obj.q0 = zeros(obj.getNumPositions, 1);
            
        end
        
        function [f, df] = dynamics(obj, t, x, u)
            
            xin = [t; x; u];
            [f,df] = dynamics_fun(obj,xin);
            fprintf('Time: %f \r', t);
            
%             df_fd = zeros(size(df));
%             step = sqrt(eps(max(xin)));
%             dxin = step*eye(length(xin));
%             for k = 1:length(xin)
%                 df_fd(:,k) = (dynamics_fun(obj, xin+dxin(:,k)) - dynamics_fun(obj, xin-dxin(:,k)))/(2*step);
%             end
% 
%             disp('Dynamics Derivative Error:');
%             disp(max(abs(df_fd(:)-df(:))));

        end
        
        function[xdot, dxdot] = dynamics_fun(obj,xin)
            
            % load in parameters
            nq = obj.getNumPositions();
            nv = obj.getNumVelocities();
            nu = obj.nu;
            nc = obj.nc;
            nl = obj.nl;
            nd = obj.nd;
            
            t = xin(1);
            q = xin(1+(1:nq));
            v = xin(1+nq+(1:nv));
            u = xin(1+nq+nv+(1:nu));
            l = xin(1+nq+nv+nu+(1:nl));
            c = xin(1+nq+nv+nu+nl+(1:nc));
            b = xin(1+nq+nv+nu+nl+nc+(1:nc*nd));
            
            % Contact basis
            kinopts = struct();
            kinopts.compute_gradients = true;
            kin = obj.doKinematics(q, v, kinopts);
            [~,~,~,~,~,~,~,~,n,D,dn,dD] = obj.contactConstraints(kin);
            
            if isempty(n)
                n = zeros(0,nq);
                dn = zeros(0,nq);
            end
            
            D = reshape(cell2mat(D')',nq,nc*nd)';
            dD = reshape(cell2mat(dD)',nq,nc*nd*nq)';
            
            % loop constraints
            K = zeros(nl, nq);
            dK = zeros(nl, nq*nq);
            
            loops = obj.loop_const;
            nli = nl/numel(loops);
            for i = 1:numel(loops)
                [~, K((i-1)*nli+(1:nli),:), dK((i-1)*nli+(1:nli),:)] = loops{i}.eval(q);
            end
            
            if isempty(dK)
                dK = zeros(0, nq);
            end
            dK = reshape(dK', nq, nl*nq)'; %
            
            % Manipulator Dynamics
            [H,C,B,dH,dC,dB] = manipulatorDynamics(obj, q, v);
            
            % form rhs and derivs
            rhs = B*u - C + K'*l + n'*c + D'*b ;
            drhs_dq = kron(u', eye(nv))*dB(:,1:nq) +  kron(c', eye(nq))*comm(nc,nq)*dn + ...
                + kron(l', eye(nv))*dK + kron(b', eye(nv))*comm(nd*nc,nv)*dD - dC(:,1:nq);
            drhs_dv = kron(u', eye(nv))*dB(:,nq+(1:nv)) - dC(:, nq+(1:nv));
            
            % form Hinv and derivs
            Hinv = inv(H);
            dH = reshape(dH,nq*nq,nq+nv);
            dHinv_dq = -kron(Hinv', Hinv)*dH(:,1:nq);
            
            % Solve for accelerations
            qddot = Hinv*rhs;
            
            % xdot and deriv
            xdot = [v; qddot];
            dxdot = [zeros(nv, 1), zeros(nv, nq), eye(nv), zeros(nv, nu+nl+nc+nc*nd); ... %dv/dt, dv/dq, dv/dv, dv/du, dv/dl, dv/dc, dv/db
                zeros(nv, 1), kron(rhs', eye(nv))*dHinv_dq + Hinv*drhs_dq, Hinv*drhs_dv, Hinv*B, Hinv*K', Hinv*n', Hinv*D']; %da/dt, da/dq, da/dv, da/du, da/dl, da/dc, da/db
            
        end
        %
        
        
        function nActuatedDOF = getNumActuatedDOF(obj)
            nActuatedDOF = numel(obj.getActuatedJoints());
        end
        
        function x0 = getInitialState(obj)
            q0 = obj.q0;
            x0 = [q0; 0*q0];
        end
    end
    
    
    
end
