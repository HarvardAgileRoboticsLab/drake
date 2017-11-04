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
            obj.P = reshape(options.p, 2, 2);
            obj.K = reshape(options.k, 2, 2);
            
            %Change gravity            
            obj = obj.setGravity(obj.grav);
            obj = compile(obj);
            
%             Add contact forces to inputs 
            [phi,~,d] = obj.contactConstraints(getZeroConfiguration(obj)); 
            nc = length(phi); obj.nc = nc;
            nd = 2*length(d); obj.nd = nd;             
            obj = obj.setNumInputs(obj.nu+nc+nd);           
                        
            % Initial state
            obj.q0 = [0.5186; 0];            
            
        end
        
        function [f, df] = dynamics(obj, t, x, u)
            
            xin = [t; x; u]; 
            [f,df] = dynamics_fun(obj,xin);
%             fprintf('Dynamics: %f \r', max(abs(f)));

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
            nd = obj.nd; 
            
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
            
            % Spring and Damper Forces
            fsd = obj.K*(q-obj.q0) + obj.P*v; 
            
            % Manipulator Dynamics
            [H,C,B,dH,dC,dB] = manipulatorDynamics(obj, q, v);
            
            % form rhs and derivs
            rhs = B*u + n'*c + D'*b - C + fsd;
            drhs_dq = kron(u', eye(nv))*dB(:,1:nq) + dn'*c + kron(b', eye(nv))*comm(nd*nc,nv)*dD - dC(:,1:nq) + obj.K;
            drhs_dv = kron(u', eye(nv))*dB(:,nq+(1:nv)) - dC(:, nq+(1:nv)) + obj.P;
            
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
