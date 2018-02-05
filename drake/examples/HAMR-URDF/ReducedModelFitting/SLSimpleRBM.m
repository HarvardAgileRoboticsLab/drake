classdef SLSimpleRBM < RigidBodyManipulator
    
    properties (SetAccess = protected, GetAccess = public)
        q0
        grav = [0; 0; -9.81e-3];
        nu              % number of actual inputs
        nc              % number of contact pairs
        nd              % number of basis vectors in polyhedral friction cone
        K               % torsional spring at universal joint (2x2 matrix)
        P               % torsional damper at universal joint (2x2 matrix)
        
    end
    
    methods
        
        function obj= SLSimpleRBM(urdf, options)
            
            typecheck(urdf,'char');
            obj = obj@RigidBodyManipulator(urdf,options);
            
            obj.nu = obj.getNumInputs();
            
            % spring damper
            for i = 1:size(options.k, 2)
                obj.K(:,:,i) = reshape(options.k(:,i), 2, 2);
            end
            for i = 1:size(options.p, 2)
                obj.P(:,:,i) = reshape(options.p(:,i), 2, 2);
            end
            
            %Change gravity
            obj = obj.setGravity(obj.grav);
            obj = compile(obj);
            
            % Add contact forces to inputs
            [phi,~,d] = obj.contactConstraints(getZeroConfiguration(obj));
            nc = length(phi); obj.nc = nc;
            nd = 2*length(d); obj.nd = nd;
            obj = obj.setNumInputs(obj.nu+nc+nd);
            
            % Initial state
            obj.q0 = [0; 0.5186];
            
        end
        
        function [f, df] = dynamics(obj, t, x, u)
            
            xin = [t; x; u];
            [f,df] = dynamics_fun(obj,xin);
            fprintf('Time: %f \r', t);
            
%             df_fd = zeros(size(df));
%             step = max(1e-8, sqrt(eps(max(abs(xin))))); 
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
            nd = obj.nd;
            nk = size(obj.K, 3);
            np = size(obj.P, 3);
            
            t = xin(1);
            q = xin(1+(1:nq));
            v = xin(1+nq+(1:nv));
            u = xin(1+nq+nv+(1:nu));
            c = xin(1+nq+nv+nu+(1:nc));
            b = xin(1+nq+nv+nu+nc+(1:nd));
            
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
            
            % Spring Forces
            fk = zeros(nq, 1);
            dfk_dq = zeros(nq);
            dq = q-obj.q0; 
            for i = 1:nk
                if i == 1
                    fk = fk + obj.K(:,:,i)*dq;
                    dfk_dq = dfk_dq + obj.K(:,:,i);
                elseif i == 2
                    fk = fk + obj.K(:,:,i)*norm(dq)^2*dq;
                    dfk_dq = dfk_dq + 3*norm(dq)^2*obj.K(:,:,i);
                else
                    error('need to define HOT correctly')
                end
            end
            
            % Damping Forces
            fp = zeros(nv, 1);
            dfp_dv = zeros(nv);
            for i = 1:np
                if i == 1
                    fp = fp + obj.P(:,:,i)*v;
                    dfp_dv = dfp_dv + obj.P(:,:,i); 
                elseif i == 2
                    fp = fp + obj.P(:,:,i)*norm(v)^2*v;
                    dfp_dv = dfp_dv + 3*norm(v)^2*obj.P(:,:,i);
                else
                    error('need to define HOT correctly')
                end
            end
            
            % Manipulator Dynamics
            [H,C,B,dH,dC,dB] = manipulatorDynamics(obj, q, v);
            
            % form rhs and derivs
            rhs = B*u + n'*c + D'*b - C + fk + fp; 
            drhs_dq = kron(u', eye(nv))*dB(:,1:nq) + dn'*c + kron(b', eye(nv))*comm(nd*nc,nv)*dD ...
                - dC(:,1:nq) + dfk_dq;
            drhs_dv = kron(u', eye(nv))*dB(:,nq+(1:nv)) - dC(:, nq+(1:nv)) + dfp_dv;
            
            % form Hinv and derivs
            Hinv = inv(H);
            dH = reshape(dH,nq*nq,nq+nv);
            dHinv_dq = -kron(Hinv', Hinv)*dH(:,1:nq);
            
            % Solve for accelerations
            qddot = Hinv*rhs;
            
            % xdot and deriv
            xdot = [v; qddot];
            dxdot = [zeros(nv, 1), zeros(nv, nq), eye(nv), zeros(nv, nu+nc+nd); ... %dv/dt, dv/dq, dv/dv, dv/du, dv/dc, dv/db
                zeros(nv, 1), kron(rhs', eye(nv))*dHinv_dq + Hinv*drhs_dq, Hinv*drhs_dv, Hinv*B, Hinv*n', Hinv*D']; %da/dt, da/dq, da/dv, da/du, da/dc, da/db
            
        end
        
        function [xdot, dxdot_dK, dxdot_dP] = dynamics_sd_fit_fun(obj, xin)
            
            
            [xdot, dxdot_dK, dxdot_dP] = dynamics_sd_fit(obj,xin);
            
            % %             load in parameters
            %             nq = obj.getNumPositions();
            %             nv = obj.getNumVelocities();
            %             nu = obj.nu;
            %
            %             dxdot_dKd = zeros(size(dxdot_dK));
            %             dxdot_dPd = zeros(size(dxdot_dP));
            %
            %             ind0 = 1+nq+nv+nu;
            %             step = sqrt(eps(max(xin(ind0:end))));
            %             dxin = step*eye(length(xin));
            %
            %             for k = ind0 +(1:nq^2)
            %                 dxdot_dKd(:, k-ind0) = (dynamics_sd_fit(obj, xin+dxin(:,k)) - ...
            %                     dynamics_sd_fit(obj, xin-dxin(:,k)))/(2*step);
            %                 dxdot_dPd(:, k-ind0) = (dynamics_sd_fit(obj, xin+dxin(:,k+nq^2)) - ...
            %                     dynamics_sd_fit(obj, xin-dxin(:,k+nq^2)))/(2*step);
            %
            %             end
            %
            %             disp('Dynamics Derivative Error:');
            %             disp(max(abs(dxdot_dKd(:)-dxdot_dK(:))));
            %             disp(max(abs(dxdot_dPd(:)-dxdot_dP(:))));
            %
        end
        
        
        function[xdot, dxdot_dK, dxdot_dP] = dynamics_sd_fit(obj,xin)
            
            % load in parameters
            nq = obj.getNumPositions();
            nv = obj.getNumVelocities();
            nu = obj.nu;
            nk = size(obj.K, 3);
            np = size(obj.P, 3);
            
            t = xin(1);
            q = xin(1+(1:nq));
            v = xin(1+nq+(1:nv));
            u = xin(1+nq+nv+(1:nu));
            K = reshape(xin(1+nq+nv+nu+(1:(nk*nq^2))), nq, nq, nk);
            P = reshape(xin(1+nq+nv+nu+nk*nq^2+(1:(np*nv^2))), nv, nv, np);
            
            % Spring Forces
            fk = zeros(nq, 1);
            dfk_dK = zeros(nk*nq, 1);
            dq = q-obj.q0;
            for i = 1:nk
                if i == 1
                    fk = fk + K(:,:,i)*dq;
                    dfk_dK(((i-1)*nq+1):(nq*i)) = dq;
                elseif i == 2
                    fk = fk + K(:,:,i)*norm(dq)^2*dq;
                    dfk_dK(((i-1)*nq+1):(nq*i)) = norm(dq)^2*dq;                   
                else
                    error('need to define HOT correctly')
                end
            end
            
            % Damping Forces
            fp = zeros(nv, 1);
            dfp_dP = zeros(np*nv, 1);
            for i = 1:np
                if i == 1
                    fp = fp + P(:,:,i)*v;
                    dfp_dP(((i-1)*nv+1):(nv*i)) = v;
                elseif i == 2
                    fp = fp + P(:,:,i)*norm(v)^2*v;
                    dfp_dP(((i-1)*nv+1):(nv*i)) = norm(v)^2*v;
                else
                    error('need to define HOT correctly')
                end
            end
            
            
            % Manipulator Dynamics
            [H,C,B,dH,~,~] = manipulatorDynamics(obj, q, v);
            
            % form rhs and derivs
            rhs = B*u - C + fk + fp;
            
            % form Hinv and derivs
            Hinv = inv(H);
            dH = reshape(dH,nq*nq,nq+nv);
            
            % Solve for accelerations
            qddot = Hinv*rhs;
            
            % xdot and deriv
            xdot = [v; qddot];
            dxdot_dK = [zeros(nq, numel(K)); kron(dfk_dK', Hinv)];
            dxdot_dP = [zeros(nq, numel(P)); kron(dfp_dP', Hinv)];
            
            
        end
        
        function obj = setK(obj, K)
            if sizecheck(K, size(obj.K))
                obj.K = K;
            else
                error('Wrongsize K')
            end
        end
        
        function obj = setP(obj, P)
            if sizecheck(P, size(obj.P))
                obj.P = P;
            else
                error('Wrongsize P')
            end
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
