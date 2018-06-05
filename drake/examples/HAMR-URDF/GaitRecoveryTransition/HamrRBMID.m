classdef HamrRBMID < HamrRBM
    
    properties (SetAccess = protected, GetAccess = public)
        
        NU              % number of actual inputs
        NC              % number of contact pairs
        ND              % number of basis vectors in polyhedral friction cone

        
    end
    
    methods
        
        function obj=HamrRBMID(urdf,options)
            
            typecheck(urdf,'char');
            
            obj = obj@HamrRBM(urdf,options);            
            obj.NU = obj.getNumInputs();            

            % Add contact forces to inputs
            [phi,~,d] = obj.contactConstraints(obj.x0(1:obj.getNumPositions()));
            NC = length(phi); obj.NC = NC;
            ND = 2*length(d); obj.ND = ND;            
           
            % Add contact and loop const as inputs
            obj = obj.setNumInputs(obj.NU + NC + NC*ND + obj.NL);
            
        end
        
        function obj = compile(obj)
            obj = compile@RigidBodyManipulator(obj);
        end
        
        
        function [f, df] = dynamics(obj, t, x, u)
            
            xin = [t; x; u];
            [f,df] = dynamics_fun(obj,xin);
%             fprintf('Time: %f \r', t);
            
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
            nu = obj.NU;
            nc = obj.NC;
            nl = obj.NL;
            nd = obj.ND;
            
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
            [~, K, dK] = obj.positionConstraints(q);
            K = K(obj.VALID_LOOPS, :);
            dK = reshape(dK(obj.VALID_LOOPS, :)', nq, nl*nq)'; %
            
            % Manipulator Dynamics
            [H,C,B,dH,dC,dB] = manipulatorDynamics(obj, q, v);
            
            % form rhs and derivs
            rhs = B*u + K'*l+ n'*c + D'*b - C;
            
            drhs_dq = kron(u', eye(nv))*dB(:,1:nq) + kron(c',eye(nq))*comm(nc,nq)*dn + ...
                kron(l', eye(nv))*dK + kron(b', eye(nv))*comm(nd*nc,nv)*dD - dC(:,1:nq);
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
            x0 = obj.x0;
        end
        
    end
end
