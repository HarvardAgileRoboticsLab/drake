classdef VariationalRigidBodyManipulator < DrakeSystem
    %This class implements a 2nd order midpoint variational integrator with
    %support for rigid body contact
    
    properties
        manip
        timestep
        twoD = false
        dirty = true
        multiple_contacts = false
        num_contact_points
        z_cache
    end
    
    methods
        function obj = VariationalRigidBodyManipulator(manipulator_or_urdf_filename,timestep,options)
            if (nargin<3)
                options=struct();
            end
            if ~isfield(options,'twoD')
                options.twoD = false;
            end
            
            typecheck(timestep,'double');
            sizecheck(timestep,1);
            
            if isempty(manipulator_or_urdf_filename) || ischar(manipulator_or_urdf_filename)
                % then make the corresponding manipulator
                w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
                if options.twoD
                    manip = PlanarRigidBodyManipulator(manipulator_or_urdf_filename,options);
                else
                    manip = RigidBodyManipulator(manipulator_or_urdf_filename,options);
                end
                warning(w);
            else
                manip = manipulator_or_urdf_filename;
            end
            
            typecheck(manip,'RigidBodyManipulator');
            obj = obj@DrakeSystem(0,manip.getNumStates(),manip.getNumInputs(),manip.getNumOutputs(),manip.isDirectFeedthrough(),manip.isTI());
            obj.manip = manip;
            if isa(manip,'PlanarRigidBodyManipulator')
                obj.twoD = true;
            end
            
            obj.z_cache = SharedDataHandle(0);
            obj.timestep = timestep;
            obj = setSampleTime(obj,[timestep;0]);
            obj = compile(obj);
            
            kin = obj.manip.doKinematics(zeros(manip.getNumPositions(),1));
            obj.num_contact_points = length(manip.contactConstraints(kin));
        end
        
        function model = compile(model)
            w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
            model.manip = model.manip.compile();
            warning(w);
            
            model = setNumDiscStates(model,model.manip.getNumContStates());
            model = setNumInputs(model,model.manip.getNumInputs());
            
            model = setInputLimits(model,model.manip.umin,model.manip.umax);
            model = setInputFrame(model,getInputFrame(model.manip));
            model = setStateFrame(model,getStateFrame(model.manip));
            model = setOutputFrame(model,getOutputFrame(model.manip));
            
            model.dirty = false;
        end
        
        function [xdn,df] = update(obj,t,x,u)
            h = obj.timestep;
            Nq = obj.manip.getNumPositions();
            Nv = obj.manip.getNumVelocities();
            Np = obj.num_contact_points;
            Nd = 4;
            
            q0 = x(1:Nq);
            v0 = x(Nq + (1:Nv));
            
            q1 = q0 + h*v0; %initial guess
            M = manipulatorDynamics(obj.manip, q0, zeros(Nv,1));
            p0 = M*v0;
            
            if Np == 0 %No contact
                r = ones(size(q0));
                while max(abs(r)) > 1e-6
                    [r,dr] = MidpointDEL(obj,p0,q0,q1);
                    dq = -dr\r;
                    alpha = 1;
                    r2 = r'*r;
                    rnew2 = r2+1;
                    while rnew2 > r2
                        q1new = q1 + alpha*dq;
                        rnew = MidpointDEL(obj,p0,q0,q1new);
                        rnew2 = rnew'*rnew;
                        alpha = alpha/2;
                    end
                    q1 = q1new;
                    r = rnew;
                end
            else %Solve with contact
                
                %z vector is stacked [q_1; c1; b1; psi; s]
                z = [q1; zeros(Np+Nd*Np+Np,1); 1];
                if length(obj.z_cache.data) == length(z)
                    z(Nq+(1:Np+Nd*Np+Np)) = obj.z_cache.data(Nq+(1:Np+Nd*Np+Np));
                end
                
                reg = 1e-3;
                [r,dr] = MidpointContact(obj,q0,p0,z,Np,Nd);
                iter = 0;
                while max(abs(r(1:(end-1)))) > 1e-6
                    iter = iter+1;
                    L = blkdiag(zeros(Nq), reg*eye(Np+Nd*Np+Np), min(100*reg,1));
                    [Q,R] = qr([dr; L],0);
                    dz = -R\(Q(1:length(r),:)'*r);
                    alpha = 1;
                    r2 = r'*r;
                    rnew2 = r2+1;
                    while rnew2 > r2 && alpha > 1e-3
                        znew = z + alpha*dz;
                        znew(Nq+(1:(Np+Nd*Np+Np))) = max(znew(Nq+(1:(Np+Nd*Np+Np))), 0);
                        znew(end) = max(znew(end), 1e-6);
                        [rnew, dr] = MidpointContact(obj,q0,p0,znew,Np,Nd);
                        rnew2 = rnew'*rnew;
                        alpha = alpha/2;
                    end
                    if alpha < .01
                        reg = min(10*reg, 1e3);
                        znew(end) = max(z(end), 1e-4);
                        [rnew, dr] = MidpointContact(obj,q0,p0,znew,Np,Nd);
                    elseif alpha > .1
                        reg = max(.5*reg, 1e-4);
                    end
                    z = znew;
                    r = rnew;
                end
                disp(iter);
                q1 = z(1:Nq);
                p1 = MidpointDLT(obj,q0,q1);
                M = manipulatorDynamics(obj.manip, q1, zeros(Nv,1));
                v1 = M\p1;
                xdn = [q1; v1];
                obj.z_cache.data = z;
            end
        end
        
        function [r,dr] = MidpointDEL(obj,p0,q0,q1)
            h = obj.timestep;
            [D1L,D2L,M] = obj.LagrangianDerivs((q0+q1)/2,(q1-q0)/h);
            r = p0 + (h/2)*D1L - D2L;
            dr = -(1/h)*M;
        end
        
        function p1 = MidpointDLT(obj,q0,q1)
            %Right Discrete Legendre transform gives momentum at end of timestep
            h = obj.timestep;
            [D1L,D2L] = obj.LagrangianDerivs((q0+q1)/2,(q1-q0)/h);
            p1 = (h/2)*D1L + D2L;
        end
        
        function [D1L,D2L,M] = LagrangianDerivs(obj,q,v)
            Nq = length(q);
            Nv = length(v);
            [M,G,~,dM] = manipulatorDynamics(obj.manip, q, zeros(Nv,1));
            dM = reshape(dM,Nq*Nq,Nq+Nv);
            dMdq = dM(:,1:Nq);
            D1L = 0.5*dMdq'*kron(v,v) - G;
            D2L = M*v;
        end
        
        function [r, dr] = MidpointContact(obj,q0,p0,z,Np,Nd)
            mu = 1; %This is currently hard coded in Drake.
            Nq = length(q0);
            h = obj.timestep;
            
            %z vector is stacked [q_1; c1; b1; psi; s]
            
            %Configurations
            q1 = z(1:Nq);
            
            %Contact force coefficients
            c = z(Nq+(1:Np));
            b = z(Nq+Np+(1:Np*Nd));
            
            %Tangential contact velocity
            psi = z(Nq+Np+Np*Nd+(1:Np));

            %Smoothing parameter
            s = z(end);
            
            %Get contact basis
            kinopts = struct();
            kinopts.compute_gradients = true;
            kin = obj.manip.doKinematics(q1, (q1-q0)/h, kinopts);
            if ~obj.manip.contact_options.use_bullet
                [phi,~,~,~,~,~,~,~,n,D,dn,dD] = obj.manip.contactConstraints(kin, obj.multiple_contacts);
                D = reshape(cell2mat(D')',Nq,Np*Nd)';
                dD = reshape(cell2mat(dD)',Nq,Np*Nd*Nq)';
            else %using bullet - this doesn't seem to work.
                kin1 = obj.manip.doKinematics(q1); %for some reason collisionDetect will crash if the kinsol has derivatives
                [phi,normal,xA,xB,idxA,idxB] = obj.manip.collisionDetect(kin1,obj.multiple_contacts);
                d = obj.manip.surfaceTangents(normal);
                [n,D,dn,dD] = obj.manip.contactConstraintDerivatives(normal,kin,idxA,idxB,xA,xB,d);
                D = reshape(cell2mat(D')',Nq,Np*Nd)';
                dD = reshape(cell2mat(dD)',Nq,Np*Nd*Nq)';
            end
            
%             %Check dn
%             delta_q = 1e-7*eye(Nq);
%             for k = 1:Nq
%                 kin = obj.manip.doKinematics(q1+delta_q(:,k));
%                 [~,~,~,~,~,~,~,~,np] = obj.manip.contactConstraints(kin, obj.multiple_contacts);
%                 kin = obj.manip.doKinematics(q1-delta_q(:,k));
%                 [~,~,~,~,~,~,~,~,nm] = obj.manip.contactConstraints(kin, obj.multiple_contacts);
%                 dn_fd(:,k) = (np(:) - nm(:))/2e-7;
%             end
%             
%             %Check dD
%             delta_q = 1e-7*eye(Nq);
%             for k = 1:Nq
%                 kin = obj.manip.doKinematics(q1+delta_q(:,k));
%                 [~,~,~,~,~,~,~,~,~,Dp] = obj.manip.contactConstraints(kin, obj.multiple_contacts);
%                 kin = obj.manip.doKinematics(q1-delta_q(:,k));
%                 [~,~,~,~,~,~,~,~,~,Dm] = obj.manip.contactConstraints(kin, obj.multiple_contacts);
%                 Dp = reshape(cell2mat(Dp')',Nq,Np*Nd)';
%                 Dm = reshape(cell2mat(Dm')',Nq,Np*Nd)';
%                 dD_fd(:,k) = (Dp(:) - Dm(:))/2e-7;
%             end
%             
%             dn = dn_fd;
%             dD = dD_fd;
            
%             dn = zeros(size(dn));
%             dD = zeros(size(dD));

            %Dynamics residual
            [r_del, dr_del] = MidpointDEL(obj,p0,q0,q1);
            r_f = r_del + h*(n'*c + D'*b);
            
            %Polyhedral friction cone
            e = ones(Nd,1);
            E = kron(eye(Np),e');
            [f1, dfa1, dfb1, dfs1] = obj.smoothFB(phi, c, s);
            [f2, dfa2, dfb2, dfs2] = obj.smoothFB(E'*phi, b, s);
            [f3, dfa3, dfb3, dfs3] = obj.smoothFB((mu*c - E*b), psi, s);
            [f4, dfa4, dfb4, dfs4] = obj.smoothFB(E'*psi + D*((q1-q0)/h), b, s);
              
            r = [r_f; f1; f2; f3; f4; exp(s)-1];
            dr = [dr_del + h*kron(c',eye(Nq))*dn + h*kron(b',eye(Nq))*dD, h*n', h*D', zeros(Nq,Np), zeros(Nq,1);
                  dfa1*n, dfb1, zeros(Np,Nd*Np+Np), dfs1;
                  dfa2*E'*n, zeros(Nd*Np,Np), dfb2, zeros(Nd*Np,Np), dfs2;
                  zeros(Np,Nq), dfa3*mu, -dfa3*E, dfb3, dfs3;
                  dfa4*(1/h)*D + kron(((q1-q0)/h)',eye(Nd*Np))*dD, zeros(Nd*Np,Np), dfb4, dfa4*E', dfs4;
                  zeros(1,Nq+Np+Nd*Np+Np), exp(s)];
        end
        
        function [f, dfda, dfdb, dfds] = smoothFB(obj,a,b,s)
            Na = length(a);
            Nb = length(b);
            if Na == Nb %vector version
                f1 = sqrt(a.*a + b.*b  + s*s*ones(Nb,1));
                f = f1 - (a + b);
                dfda = diag(a./f1) - eye(Na);
                dfdb = diag(b./f1) - eye(Nb);
                dfds = s*ones(Nb,1)./f1;
            elseif Na == 1
                f1 = sqrt(a*a*ones(Nb,1) + b.*b  + s*s*ones(Nb,1));
                f = f1 - (a*ones(Nb,1) + b);
                dfda = a*ones(Nb,1)./f1 - ones(Nb,1);
                dfdb = diag(b./f1) - eye(Nb);
                dfds = s*ones(Nb,1)./f1;
            else %Nb == 1
                f1 = sqrt(a.*a + b*b*ones(Na,1)  + s*s*ones(Na,1));
                f = f1 - (a + b*ones(Na,1));
                dfda = diag(a./f1) - ones(Na,1);
                dfdb = b*ones(Na,1)./f1 - eye(Nb);
                dfds = s*ones(Na,1)./f1;
            end
        end
        
        function x0 = getInitialState(obj)
            if ~isempty(obj.initial_state)
                x0 = obj.initial_state;
                return;
            end
            
            x0 = obj.manip.getInitialState();
%             for i=1:length(obj.sensor)
%                 if isa(obj.sensor{i},'TimeSteppingRigidBodySensorWithState')
%                     x0 = [x0; obj.sensor{i}.getInitialState(obj)];
%                 end
%             end
        end
        
        function y = output(obj,t,x,u)
            y = obj.manip.output(t,x,u);
        end
        
    end
    
    
    methods  % pass through methods (to the manipulator)
        function B = getB(obj)
            B = getB(obj.manip);
        end
        
        function g = getGravity(obj)
            g = getGravity(obj.manip);
        end
        
        function num_q = getNumPositions(obj)
            num_q = obj.manip.num_positions;
        end
        
        function num_v = getNumVelocities(obj)
            num_v = obj.manip.getNumVelocities();
        end
        
        function obj = setStateFrame(obj,fr)
            obj = setStateFrame@DrakeSystem(obj,fr);
            
            % make sure there is a transform defined to and from the
            % manipulator state frame.  (the trivial transform is the correct
            % one)
            if ~isempty(obj.manip) % this also gets called on the initial constructor
                mfr = getStateFrame(obj.manip);
                if isempty(findTransform(fr,mfr))
                    addTransform(fr,AffineTransform(fr,mfr,eye(obj.manip.num_x,obj.num_x),zeros(obj.manip.num_x,1)));
                end
                if isempty(findTransform(mfr,fr))
                    addTransform(mfr,AffineTransform(mfr,fr,eye(obj.num_x,obj.manip.num_x),zeros(obj.num_x,1)));
                end
            end
        end
        
        function obj = setTerrain(obj,varargin)
            obj.manip = setTerrain(obj.manip,varargin{:});
        end
        
        function terrain = getTerrain(obj)
            terrain = obj.manip.terrain;
        end
        
        function varargout = getTerrainHeight(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = getTerrainHeight(obj.manip,varargin{:});
        end
        
        function obj = setJointLimits(obj,varargin)
            obj.manip = setJointLimits(obj.manip,varargin{:});
        end
        
        function obj=addRobotFromURDF(obj,varargin)
            if obj.twoD
                w = warning('off','Drake:PlanarRigidBodyManipulator:UnsupportedContactPoints');
                warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
            else
                w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
            end
            obj.manip=addRobotFromURDF(obj.manip,varargin{:});
            obj=compile(obj);  % note: compiles the manip twice, but it's ok.
            warning(w);
        end
        
        function obj=addRobotFromSDF(obj,varargin)
            obj.manip=addRobotFromSDF(obj.manip,varargin{:});
            obj=compile(obj);  % note: compiles the manip twice, but it's ok.
        end
        
        function varargout = doKinematics(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}]=doKinematics(obj.manip,varargin{:});
        end
        
        function varargout = forwardKin(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}]=forwardKin(obj.manip,varargin{:});
        end
        
        function varargout = bodyKin(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}]=bodyKin(obj.manip,varargin{:});
        end
        
        function varargout = approximateIK(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}]=approximateIK(obj.manip,varargin{:});
        end
        
        function varargout = inverseKin(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}]=inverseKin(obj.manip,varargin{:});
        end
        
        function varargout = inverseKinPointwise(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = inverseKinPointwise(obj.manip,varargin{:});
        end
        
        function varargout = inverseKinTraj(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = inverseKinTraj(obj.manip,varargin{:});
        end
        
        function varargout = inverseKinWrapup(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = inverseKinWrapup(obj.manip,varargin{:});
        end
        
        function varargout = findFixedPoint(obj,x0,varargin)
            varargout=cell(1,nargout);
            if isnumeric(x0)
                x0 = Point(obj.getStateFrame(),x0);
            end
            [varargout{:}]=findFixedPoint(obj.manip,x0,varargin{:});
            varargout{1} = varargout{1}.inFrame(obj.getStateFrame());
        end
        
        function varargout = collisionDetect(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}]=collisionDetect(obj.manip,varargin{:});
        end
        
        function varargout = collisionDetectTerrain(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}]=collisionDetectTerrain(obj.manip,varargin{:});
        end
        
        function [obj,id] = addStateConstraint(obj,con)
            % keep two copies of the constraints around ... :(
            % todo: re-evaluate whether that is really necessary
            [obj,id] = addStateConstraint@DrakeSystem(obj,con);
            [obj.manip,manip_id] = obj.manip.addStateConstraint(obj,con);
            assert(id==manip_id);
        end
        
        function obj = updateStateConstraint(obj,id,con)
            obj = updateStateConstraint@DrakeSystem(obj,id,con);
            obj.manip = updateStateConstraint(obj.manip,id,con);
        end
        
        function obj = removeAllStateConstraints(obj)
            obj = removeAllStateConstraints@DrakeSystem(obj);
            obj.manip = removeAllStateConstraints(obj.manip);
        end
        
        function varargout = positionConstraints(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = positionConstraints(obj.manip,varargin{:});
        end
        
        function varargout = velocityConstraints(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = velocityConstraints(obj.manip,varargin{:});
        end
        
        function varargout = manipulatorDynamics(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = manipulatorDynamics(obj.manip,varargin{:});
        end
        
        function varargout = contactConstraints(obj,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = contactConstraints(obj.manip,varargin{:});
        end
        
        function varargout = contactConstraintsBV(obj,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = contactConstraintsBV(obj.manip,varargin{:});
        end
        
        function varargout = pairwiseContactConstraints(obj,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = pairwiseContactConstraints(obj.manip,varargin{:});
        end
        
        function varargout = pairwiseContactConstraintsBV(obj,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = pairwiseContactConstraintsBV(obj.manip,varargin{:});
        end
        
        function varargout = resolveConstraints(obj,x0,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = resolveConstraints(obj.manip,x0,varargin{:});
            varargout{1} = varargout{1}.inFrame(obj.getStateFrame());
        end
        
        function varargout = getMass(obj,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = getMass(obj.manip,varargin{:});
        end
        
        function varargout = getCOM(obj,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = getCOM(obj.manip,varargin{:});
        end
        
        function varargout = centerOfMassJacobianDotTimesV(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = centerOfMassJacobianDotTimesV(obj.manip,varargin{:});
        end
        
        function varargout = centroidalMomentumMatrixDotTimesV(obj,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = centroidalMomentumMatrixDotTimesV(obj.manip,varargin{:});
        end
        
        function varargout = centroidalMomentumMatrix(obj,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = centroidalMomentumMatrix(obj.manip,varargin{:});
        end
        
        function varargout = parseBodyOrFrameID(obj,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = parseBodyOrFrameID(obj.manip,varargin{:});
        end
        
        function joint_ind = findJointId(model,varargin)
            joint_ind = findJointId(model.manip,varargin{:});
        end
        
        function body_ind = findLinkId(model,varargin)
            body_ind = findLinkId(model.manip,varargin{:});
        end
        
        function indices = findPositionIndices(model, varargin)
            indices = findPositionIndices(model.manip,varargin{:});
        end
        
        function body = findLink(model,varargin)
            body = findLink(model.manip,varargin{:});
        end
        
        function frame_id = findFrameId(model,varargin)
            frame_id = findFrameId(model.manip,varargin{:});
        end
        
        function ancestor_bodies = findAncestorBodies(obj, body_index)
            ancestor_bodies = obj.manip.findAncestorBodies(body_index);
        end
        
        function [body_path, joint_path, signs] = findKinematicPath(obj, start_body, end_body)
            [body_path, joint_path, signs] = obj.manip.findKinematicPath(start_body, end_body);
        end
        
        function obj = weldJoint(obj,body_ind_or_joint_name,robot)
            if nargin>2
                obj.manip = weldJoint(obj.manip,body_ind_or_joint_name,robot);
            else
                obj.manip = weldJoint(obj.manip,body_ind_or_joint_name);
            end
            obj.dirty = true;
        end
        
        function body = getBody(model,varargin)
            body = getBody(model.manip,varargin{:});
        end
        
        function frame = getFrame(model,varargin)
            frame = getFrame(model.manip,varargin{:});
        end
        
        function str = getBodyOrFrameName(obj,varargin)
            str = obj.manip.getBodyOrFrameName(varargin{:});
        end
        
        function model = setBody(model,varargin)
            model.manip = setBody(model.manip,varargin{:});
            model.dirty = true;
        end
        
        function v = constructVisualizer(obj,varargin)
            v = constructVisualizer(obj.manip,varargin{:});
        end
        
        function getNumContacts(~)
            error('getNumContacts is no longer supported, in anticipation of alowing multiple contacts per body pair. Use getNumContactPairs for cases where the number of contacts is fixed');
        end
        
        function n=getNumContactPairs(obj)
            n = obj.manip.getNumContactPairs;
        end
        
        function c = getBodyContacts(obj,body_idx)
            c = obj.manip.body(body_idx).collision_geometry;
        end
        
        function addContactShapeToBody(varargin)
            errorDeprecatedFunction('addCollisionGeometryToBody');
        end
        
        function obj = addCollisionGeometryToBody(obj,varargin)
            obj.manip = addCollisionGeometryToBody(obj.manip,varargin{:});
        end
        
        function addVisualShapeToBody(varargin)
            errorDeprecatedFunction('addVisualGeometryToBody');
        end
        
        function obj = addVisualGeometryToBody(obj,varargin)
            obj.manip = addVisualGeometryToBody(obj.manip,varargin{:});
        end
        
        function addShapeToBody(varargin)
            errorDeprecatedFunction('addGeometryToBody');
        end
        
        function obj = addGeometryToBody(obj,varargin)
            obj.manip = addGeometryToBody(obj.manip,varargin{:});
        end
        
        function replaceContactShapesWithCHull(varargin)
            errorDeprecatedFunction('replaceCollisionGeometryWithConvexHull');
        end
        
        function obj = replaceCollisionGeometryWithConvexHull(obj,body_indices,varargin)
            obj.manip = replaceCollisionGeometryWithConvexHull(obj.manip,body_indices,varargin{:});
        end
        
        function getContactShapeGroupNames(varargin)
            errorDeprecatedFunction('getCollisionGeometryGroupNames');
        end
        
        function groups = getCollisionGeometryGroupNames(obj)
            groups = getCollisionGeometryGroupNames(obj.manip);
        end
        
        function f_friction = computeFrictionForce(obj,qd)
            f_friction = computeFrictionForce(obj.manip,qd);
        end
        
        function obj = removeCollisionGroups(obj,contact_groups)
            obj.manip = removeCollisionGroups(obj.manip,contact_groups);
        end
        
        function obj = removeCollisionGroupsExcept(obj,varargin)
            obj.manip = removeCollisionGroupsExcept(obj.manip,varargin{:});
        end
        
        function str = getLinkName(obj,body_ind)
            str = obj.manip.getLinkName(body_ind);
        end
        
        function link_names = getLinkNames(obj)
            link_names =  {obj.manip.body.linkname}';
        end
        
        function joint_names = getJointNames(obj)
            joint_names =  {obj.manip.body.jointname}';
        end
        
        function num_bodies = getNumBodies(obj)
            num_bodies = length(obj.manip.body);
        end
        
        function [jl_min, jl_max] = getJointLimits(obj)
            jl_min = obj.manip.joint_limit_min;
            jl_max = obj.manip.joint_limit_max;
        end
        
        function varargout = jointLimitConstraints(obj,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = jointLimitConstraints(obj.manip,varargin{:});
        end
        
        function index = getActuatedJoints(obj)
            index = getActuatedJoints(obj.manip);
        end
        
        function ptr = getMexModelPtr(obj)
            ptr = getMexModelPtr(obj.manip);
        end
        
        function [phi,Jphi] = closestDistance(obj,varargin)
            [phi,Jphi] = closestDistance(obj.manip,varargin{:});
        end
        
        function obj = addLinksToCollisionFilterGroup(obj,linknames,collision_fg_name,robotnums)
            obj.manip = addLinksToCollisionFilterGroup(obj.manip,linknames,collision_fg_name,robotnums);
        end
        
        function out = name(obj)
            out = obj.manip.name;
        end
        
        function fr = getParamFrame(model)
            fr = getParamFrame(model.manip);
        end
        
        function model = setParams(model,p)
            model.manip = setParams(model.manip,p);
        end
        
        function terrain_contact_point_struct = getTerrainContactPoints(obj,varargin)
            terrain_contact_point_struct = getTerrainContactPoints(obj.manip,varargin{:});
        end
        
        function varargout = terrainContactPositions(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = terrainContactPositions(obj.manip,varargin{:});
        end
        
        function varargout = terrainContactJacobianDotTimesV(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = terrainContactJacobianDotTimesV(obj.manip,varargin{:});
        end
        
        function distance = collisionRaycast(obj, kinsol, origin, point_on_ray, use_margins)
            if nargin < 5
                use_margins = true;
            end
            distance = collisionRaycast(obj.manip, kinsol, origin, point_on_ray, use_margins);
        end
        
        
    end
    
end

