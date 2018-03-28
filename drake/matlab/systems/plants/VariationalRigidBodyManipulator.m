classdef VariationalRigidBodyManipulator < DrakeSystem
    %This class implements a 2nd order midpoint variational integrator with
    %support for rigid body contact
    
    properties
        plant
        timestep
        twoD = false
        dirty = true
        multiple_contacts = false
        nQ
        nV
        nU
        nD
        nC
        angle_inds
        options
    end
    
    methods
        function obj = VariationalRigidBodyManipulator(manipulator_or_urdf_filename,timestep,options)
            if (nargin<3)
                options=struct();
            end
            if ~isfield(options,'active_collision_options')
                options.active_collision_options.terrain_only = true;
            end
            if ~isfield(options,'integration_method')
                options.integration_method = VariationalTrajectoryOptimization.MIDPOINT;
            end
            if ~isfield(options,'multiple_contacts')
                options.multiple_contacts = false;
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
                    plant = PlanarRigidBodyManipulator(manipulator_or_urdf_filename,options);
                else
                    plant = RigidBodyManipulator(manipulator_or_urdf_filename,options);
                end
                warning(w);
            else
                plant = manipulator_or_urdf_filename;
            end
            
            typecheck(plant,'RigidBodyManipulator');
            obj = obj@DrakeSystem(0,plant.getNumStates(),plant.getNumInputs(),plant.getNumOutputs(),plant.isDirectFeedthrough(),plant.isTI());
            obj.plant = plant;
            if isa(plant,'PlanarRigidBodyManipulator')
                obj.twoD = true;
            end
            
            %Find indices that correspond to joint angles so we can
            %watch out for wrap-around issues
            obj.angle_inds = [obj.plant.body(~[obj.plant.body.floating]).position_num]';
            obj.angle_inds = obj.angle_inds(obj.angle_inds>0);
            
            obj.timestep = timestep;
            obj = setSampleTime(obj,[timestep;0]);
            obj = compile(obj);
            
            obj.nQ = plant.getNumPositions();
            obj.nV = plant.getNumVelocities();
            obj.nU = plant.getNumInputs();
            kin = obj.plant.doKinematics(zeros(plant.getNumPositions(),1));
            obj.nC = length(plant.contactConstraints(kin));
            obj.nD = 4;
            obj.options = options;
            
        end
        
        function model = compile(model)
            w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
            model.plant = model.plant.compile();
            warning(w);
            
            model = setNumDiscStates(model,model.plant.getNumContStates());
            model = setNumInputs(model,model.plant.getNumInputs());
            
            model = setInputLimits(model,model.plant.umin,model.plant.umax);
            model = setInputFrame(model,getInputFrame(model.plant));
            model = setStateFrame(model,getStateFrame(model.plant));
            model = setOutputFrame(model,getOutputFrame(model.plant));
            
            model.dirty = false;
        end
        
        function [xn,dxn] = update(obj,t,x,u)
            h = obj.timestep;
            nD = obj.nD;
            nC = obj.nC;
            nQ = obj.nQ;
            nV = obj.nV;
            nU = obj.nU;
            
            q0 = x(1:nQ);
            v0 = x(nQ + (1:nV));
            
            M = manipulatorDynamics(obj.plant, q0, zeros(nV,1));
            p0 = M*v0;
            
            %initial guess
            xg = [q0; zeros(nC,1); zeros(nC*nD,1); zeros(nC,1); zeros(nC*nD,1); 1];
            
            %Bounds
            lb = [-inf*ones(nQ,1); zeros(2*nC*(1+nD),1); 1e-6];
            ub = inf*ones(1+nQ+2*nC*(1+nD),1);
            
            opts = optimoptions('fmincon','Display','none','Algorithm','sqp','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'StepTolerance',1e-12);%,'FiniteDifferenceType','central','CheckGradients',true);
            xopt = fmincon(@(x)MidpointObj(obj,x),xg,[],[],[],[],lb,ub,@(x)MidpointConstr(obj,q0,p0,u,x),opts);
            
            q1 = xopt(1:nQ);
            c = xopt(nQ+(1:nC));
            b = xopt(nQ+nC+(1:(nD*nC)));
            qm = obj.qavg(q0,q1);
            vm = obj.qdiff(q0,q1,h);
            [D1L,D2L,~,~,~,B,~] = obj.LagrangianDerivs(qm,vm);
            %Contact basis
            kinopts = struct();
            kinopts.compute_gradients = true;
            kin = obj.plant.doKinematics(q1, vm, kinopts);
            [~,~,~,~,~,~,~,~,n,D] = obj.plant.contactConstraints(kin, obj.options.multiple_contacts, obj.options.active_collision_options);
            if isempty(n)
                n = zeros(0,nQ);
            end
            D = reshape(cell2mat(D')',nQ,nC*nD)';

            p1 = D2L + (h/2)*D1L + h*(B*u + n'*c + D'*b);
            M = manipulatorDynamics(obj.plant, q1, zeros(nV,1));
            v1 = M\p1;
            xn = [q1; v1];
        end
        
        function [f,df] = MidpointObj(obj,x)
            s = x(end);
            f = exp(s)-1;
            df = [zeros(1,length(x)-1), exp(s)];
        end
        
        function [cin,ceq,dcin,dceq] = MidpointConstr(obj,q0,p0,u,x)
            h = obj.timestep;
            nD = obj.nD;
            nC = obj.nC;
            nQ = obj.nQ;
            nV = obj.nV;
            nU = obj.nU;
            
            mu = 1;
            
            q1 = x(1:nQ);
            c = x(nQ+(1:nC));
            b = x(nQ+nC+(1:(nD*nC)));
            psi = x(nQ+nC+(nD*nC)+(1:nC));
            eta = x(nQ+nC+(nD*nC)+nC+(1:(nD*nC)));
            s = x(end);
            
            qm = obj.qavg(q0,q1);
            dqm = 0.5*eye(nQ);
            vm = obj.qdiff(q0,q1,h);
            dvm = (1/h)*eye(nQ);
            
            [D1L,D2L,D1D1L,D1D2L,D2D2L,B,dB] = obj.LagrangianDerivs(qm,vm);
            
            %Contact basis
            kinopts = struct();
            kinopts.compute_gradients = true;
            kin = obj.plant.doKinematics(q1, vm, kinopts);
            [phi,~,~,~,~,~,~,~,n,D,dn,dD] = obj.plant.contactConstraints(kin, obj.options.multiple_contacts, obj.options.active_collision_options);
            if isempty(n)
                n = zeros(0,nQ);
                dn = zeros(0,nQ);
            end
            D = reshape(cell2mat(D')',nQ,nC*nD)';
            dD = reshape(cell2mat(dD)',nQ,nC*nD*nQ)';
            
            e = ones(nD,1);
            E = kron(eye(nC),e');
            
            %Discrete Euler-Lagrange equation
            ceq1 = -p0 - (h/2)*D1L + D2L - h*(B*u + n'*c + D'*b); % = 0
            dceq1 = [(-h/2)*(D1D1L*dqm + D1D2L'*dvm) + D1D2L*dqm + D2D2L*dvm + -h*kron(u',eye(nQ))*dB - h*kron(c',eye(nQ))*comm(nC,nQ)*dn - h*kron(b',eye(nQ))*comm(nC*nD,nQ)*dD, -h*n', -h*D', zeros(nQ,1+nC*(1+nD))];
           
            %Tangential velocity
            ceq2 = D*vm + E'*psi - eta; % = 0
            dceq2 = [D*dvm + kron(vm',eye(nC*nD))*dD, zeros(nC*nD, nC*(1+nD)), E', -eye(nC*nD), zeros(nC*nD,1)];
            
            ceq = [ceq1; ceq2];
            dceq = [dceq1; dceq2]';
            
            %Signed distance
            cin1 = -phi; % <= 0
            dcin1 = [-n, zeros(nC,2*nC*(1+nD)+1)];
            
            %Friction cone
            cin2 = E*b - mu*c; % <= 0
            dcin2 = [zeros(nC, nQ), -mu*eye(nC), E, zeros(nC,nC*(1+nD)+1)];
            
            %Normal force complementarity
            cin3 = phi'*c - s; % <= 0
            dcin3 = [c'*n, phi', zeros(1,nC*nD+nC+nC*nD), -1];
                
            %Tangential velocity complementarity
            cin4 = (mu*c - E*b)'*psi - s; % <= 0
            dcin4 = [zeros(1,nQ), mu*psi', -psi'*E, (mu*c - E*b)', zeros(1,nD*nC), -1];
            
            %Friction complementarity
            cin5 = eta'*b - s; % <= 0
            dcin5 = [zeros(1,nQ), zeros(1,nC), eta', zeros(1,nC), b', -1];
            
            cin = [cin1; cin2; cin3; cin4; cin5];
            dcin = [dcin1; dcin2; dcin3; dcin4; dcin5]';
            
        end
        
        
        function [r,drdq1] = MidpointDEL(obj,p0,q0,q1)
            h = obj.timestep;
            [D1L,D2L,D1D1L,D1D2L,D2D2L] = obj.LagrangianDerivs((q0+q1)/2,(q1-q0)/h);
            r = p0 + (h/2)*D1L - D2L;
            
            drdq1 = (h/2)*((1/2)*D1D1L + (1/h)*D1D2L') - ((1/2)*D1D2L + (1/h)*D2D2L);
        end
        
        function p1 = MidpointRightDLT(obj,q0,q1)
            %Right Discrete Legendre transform gives momentum at end of timestep
            h = obj.timestep;
            [D1L,D2L] = obj.LagrangianDerivs(obj.qavg(q0,q1),obj.qdiff(q0,q1,h));
            p1 = (h/2)*D1L + D2L;
        end
        
        function p0 = MidpointLeftDLT(obj,q0,q1)
            %Left Discrete Legendre transform gives momentum at end of timestep
            h = obj.timestep;
            [D1L,D2L] = obj.LagrangianDerivs(obj.qavg(q0,q1),obj.qdiff(q0,q1,h));
            p0 = -(h/2)*D1L + D2L;
        end
        
        function [D1L,D2L,D1D1L,D1D2L,D2D2L,B,dBdq] = LagrangianDerivs(obj,q2,v)
            nq = length(q2);
            nv = length(v);
            
            [M,G,B,dM,dG,dB] = obj.plant.manipulatorDynamics(q2, zeros(nv,1));
            dM = reshape(dM,nq*nq,nq+nv);
            dMdq = dM(:,1:nq);

            [~,d2T] = obj.plant.kineticEnergyDerivatives(q2, v);
            
            dBdq = dB(:,1:nq);
            
            D1L = 0.5*dMdq'*kron(v,v) - G; 
            D2L = M*v;
            
            D1D1L = d2T(1:nq,1:nq) - dG(:,1:nq);
            D1D2L = kron(v',eye(nq))*dMdq;
            D2D2L = M;
        end
        

        
        function qm = qavg(obj,q1,q2)
            qm = (q1+q2)/2;
            if ~isempty(obj.angle_inds)
                qm(obj.angle_inds) = angleAverage(q1(obj.angle_inds),q2(obj.angle_inds));
            end
        end
        
        function vm = qdiff(obj,q1,q2,h)
            vm = (q2-q1)/h;
            vm(obj.angle_inds) = angleDiff(q1(obj.angle_inds),q2(obj.angle_inds))/h;
        end
        
        function x0 = getInitialState(obj)
            if ~isempty(obj.initial_state)
                x0 = obj.initial_state;
                return;
            end
            
            x0 = obj.plant.getInitialState();
%             for i=1:length(obj.sensor)
%                 if isa(obj.sensor{i},'TimeSteppingRigidBodySensorWithState')
%                     x0 = [x0; obj.sensor{i}.getInitialState(obj)];
%                 end
%             end
        end
        
        function y = output(obj,t,x,u)
            y = obj.plant.output(t,x,u);
        end
        
    end
    
    
    methods  % pass through methods (to the manipulator)
        function B = getB(obj)
            B = getB(obj.plant);
        end
        
        function g = getGravity(obj)
            g = getGravity(obj.plant);
        end
        
        function num_q = getNumPositions(obj)
            num_q = obj.plant.num_positions;
        end
        
        function num_v = getNumVelocities(obj)
            num_v = obj.plant.getNumVelocities();
        end
        
        function obj = setStateFrame(obj,fr)
            obj = setStateFrame@DrakeSystem(obj,fr);
            
            % make sure there is a transform defined to and from the
            % manipulator state frame.  (the trivial transform is the correct
            % one)
            if ~isempty(obj.plant) % this also gets called on the initial constructor
                mfr = getStateFrame(obj.plant);
                if isempty(findTransform(fr,mfr))
                    addTransform(fr,AffineTransform(fr,mfr,eye(obj.plant.num_x,obj.num_x),zeros(obj.plant.num_x,1)));
                end
                if isempty(findTransform(mfr,fr))
                    addTransform(mfr,AffineTransform(mfr,fr,eye(obj.num_x,obj.plant.num_x),zeros(obj.num_x,1)));
                end
            end
        end
        
        function obj = setTerrain(obj,varargin)
            obj.plant = setTerrain(obj.plant,varargin{:});
        end
        
        function terrain = getTerrain(obj)
            terrain = obj.plant.terrain;
        end
        
        function varargout = getTerrainHeight(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = getTerrainHeight(obj.plant,varargin{:});
        end
        
        function obj = setJointLimits(obj,varargin)
            obj.plant = setJointLimits(obj.plant,varargin{:});
        end
        
        function obj=addRobotFromURDF(obj,varargin)
            if obj.twoD
                w = warning('off','Drake:PlanarRigidBodyManipulator:UnsupportedContactPoints');
                warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
            else
                w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
            end
            obj.plant=addRobotFromURDF(obj.plant,varargin{:});
            obj=compile(obj);  % note: compiles the plant twice, but it's ok.
            warning(w);
        end
        
        function obj=addRobotFromSDF(obj,varargin)
            obj.plant=addRobotFromSDF(obj.plant,varargin{:});
            obj=compile(obj);  % note: compiles the plant twice, but it's ok.
        end
        
        function varargout = doKinematics(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}]=doKinematics(obj.plant,varargin{:});
        end
        
        function varargout = forwardKin(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}]=forwardKin(obj.plant,varargin{:});
        end
        
        function varargout = bodyKin(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}]=bodyKin(obj.plant,varargin{:});
        end
        
        function varargout = approximateIK(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}]=approximateIK(obj.plant,varargin{:});
        end
        
        function varargout = inverseKin(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}]=inverseKin(obj.plant,varargin{:});
        end
        
        function varargout = inverseKinPointwise(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = inverseKinPointwise(obj.plant,varargin{:});
        end
        
        function varargout = inverseKinTraj(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = inverseKinTraj(obj.plant,varargin{:});
        end
        
        function varargout = inverseKinWrapup(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = inverseKinWrapup(obj.plant,varargin{:});
        end
        
        function varargout = findFixedPoint(obj,x0,varargin)
            varargout=cell(1,nargout);
            if isnumeric(x0)
                x0 = Point(obj.getStateFrame(),x0);
            end
            [varargout{:}]=findFixedPoint(obj.plant,x0,varargin{:});
            varargout{1} = varargout{1}.inFrame(obj.getStateFrame());
        end
        
        function varargout = collisionDetect(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}]=collisionDetect(obj.plant,varargin{:});
        end
        
        function varargout = collisionDetectTerrain(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}]=collisionDetectTerrain(obj.plant,varargin{:});
        end
        
        function [obj,id] = addStateConstraint(obj,con)
            % keep two copies of the constraints around ... :(
            % todo: re-evaluate whether that is really necessary
            [obj,id] = addStateConstraint@DrakeSystem(obj,con);
            [obj.plant,manip_id] = obj.plant.addStateConstraint(obj,con);
            assert(id==manip_id);
        end
        
        function obj = updateStateConstraint(obj,id,con)
            obj = updateStateConstraint@DrakeSystem(obj,id,con);
            obj.plant = updateStateConstraint(obj.plant,id,con);
        end
        
        function obj = removeAllStateConstraints(obj)
            obj = removeAllStateConstraints@DrakeSystem(obj);
            obj.plant = removeAllStateConstraints(obj.plant);
        end
        
        function varargout = positionConstraints(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = positionConstraints(obj.plant,varargin{:});
        end
        
        function varargout = velocityConstraints(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = velocityConstraints(obj.plant,varargin{:});
        end
        
        function varargout = manipulatorDynamics(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = manipulatorDynamics(obj.plant,varargin{:});
        end
        
        function varargout = contactConstraints(obj,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = contactConstraints(obj.plant,varargin{:});
        end
        
        function varargout = contactConstraintsBV(obj,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = contactConstraintsBV(obj.plant,varargin{:});
        end
        
        function varargout = pairwiseContactConstraints(obj,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = pairwiseContactConstraints(obj.plant,varargin{:});
        end
        
        function varargout = pairwiseContactConstraintsBV(obj,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = pairwiseContactConstraintsBV(obj.plant,varargin{:});
        end
        
        function varargout = resolveConstraints(obj,x0,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = resolveConstraints(obj.plant,x0,varargin{:});
            varargout{1} = varargout{1}.inFrame(obj.getStateFrame());
        end
        
        function varargout = getMass(obj,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = getMass(obj.plant,varargin{:});
        end
        
        function varargout = getCOM(obj,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = getCOM(obj.plant,varargin{:});
        end
        
        function varargout = centerOfMassJacobianDotTimesV(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = centerOfMassJacobianDotTimesV(obj.plant,varargin{:});
        end
        
        function varargout = centroidalMomentumMatrixDotTimesV(obj,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = centroidalMomentumMatrixDotTimesV(obj.plant,varargin{:});
        end
        
        function varargout = centroidalMomentumMatrix(obj,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = centroidalMomentumMatrix(obj.plant,varargin{:});
        end
        
        function varargout = parseBodyOrFrameID(obj,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = parseBodyOrFrameID(obj.plant,varargin{:});
        end
        
        function joint_ind = findJointId(model,varargin)
            joint_ind = findJointId(model.plant,varargin{:});
        end
        
        function body_ind = findLinkId(model,varargin)
            body_ind = findLinkId(model.plant,varargin{:});
        end
        
        function indices = findPositionIndices(model, varargin)
            indices = findPositionIndices(model.plant,varargin{:});
        end
        
        function body = findLink(model,varargin)
            body = findLink(model.plant,varargin{:});
        end
        
        function frame_id = findFrameId(model,varargin)
            frame_id = findFrameId(model.plant,varargin{:});
        end
        
        function ancestor_bodies = findAncestorBodies(obj, body_index)
            ancestor_bodies = obj.plant.findAncestorBodies(body_index);
        end
        
        function [body_path, joint_path, signs] = findKinematicPath(obj, start_body, end_body)
            [body_path, joint_path, signs] = obj.plant.findKinematicPath(start_body, end_body);
        end
        
        function obj = weldJoint(obj,body_ind_or_joint_name,robot)
            if nargin>2
                obj.plant = weldJoint(obj.plant,body_ind_or_joint_name,robot);
            else
                obj.plant = weldJoint(obj.plant,body_ind_or_joint_name);
            end
            obj.dirty = true;
        end
        
        function body = getBody(model,varargin)
            body = getBody(model.plant,varargin{:});
        end
        
        function frame = getFrame(model,varargin)
            frame = getFrame(model.plant,varargin{:});
        end
        
        function str = getBodyOrFrameName(obj,varargin)
            str = obj.plant.getBodyOrFrameName(varargin{:});
        end
        
        function model = setBody(model,varargin)
            model.plant = setBody(model.plant,varargin{:});
            model.dirty = true;
        end
        
        function v = constructVisualizer(obj,varargin)
            v = constructVisualizer(obj.plant,varargin{:});
        end
        
        function getNumContacts(~)
            error('getNumContacts is no longer supported, in anticipation of alowing multiple contacts per body pair. Use getNumContactPairs for cases where the number of contacts is fixed');
        end
        
        function n=getNumContactPairs(obj)
            n = obj.plant.getNumContactPairs;
        end
        
        function c = getBodyContacts(obj,body_idx)
            c = obj.plant.body(body_idx).collision_geometry;
        end
        
        function addContactShapeToBody(varargin)
            errorDeprecatedFunction('addCollisionGeometryToBody');
        end
        
        function obj = addCollisionGeometryToBody(obj,varargin)
            obj.plant = addCollisionGeometryToBody(obj.plant,varargin{:});
        end
        
        function addVisualShapeToBody(varargin)
            errorDeprecatedFunction('addVisualGeometryToBody');
        end
        
        function obj = addVisualGeometryToBody(obj,varargin)
            obj.plant = addVisualGeometryToBody(obj.plant,varargin{:});
        end
        
        function addShapeToBody(varargin)
            errorDeprecatedFunction('addGeometryToBody');
        end
        
        function obj = addGeometryToBody(obj,varargin)
            obj.plant = addGeometryToBody(obj.plant,varargin{:});
        end
        
        function replaceContactShapesWithCHull(varargin)
            errorDeprecatedFunction('replaceCollisionGeometryWithConvexHull');
        end
        
        function obj = replaceCollisionGeometryWithConvexHull(obj,body_indices,varargin)
            obj.plant = replaceCollisionGeometryWithConvexHull(obj.plant,body_indices,varargin{:});
        end
        
        function getContactShapeGroupNames(varargin)
            errorDeprecatedFunction('getCollisionGeometryGroupNames');
        end
        
        function groups = getCollisionGeometryGroupNames(obj)
            groups = getCollisionGeometryGroupNames(obj.plant);
        end
        
        function f_friction = computeFrictionForce(obj,qd)
            f_friction = computeFrictionForce(obj.plant,qd);
        end
        
        function obj = removeCollisionGroups(obj,contact_groups)
            obj.plant = removeCollisionGroups(obj.plant,contact_groups);
        end
        
        function obj = removeCollisionGroupsExcept(obj,varargin)
            obj.plant = removeCollisionGroupsExcept(obj.plant,varargin{:});
        end
        
        function str = getLinkName(obj,body_ind)
            str = obj.plant.getLinkName(body_ind);
        end
        
        function link_names = getLinkNames(obj)
            link_names =  {obj.plant.body.linkname}';
        end
        
        function joint_names = getJointNames(obj)
            joint_names =  {obj.plant.body.jointname}';
        end
        
        function num_bodies = getNumBodies(obj)
            num_bodies = length(obj.plant.body);
        end
        
        function [jl_min, jl_max] = getJointLimits(obj)
            jl_min = obj.plant.joint_limit_min;
            jl_max = obj.plant.joint_limit_max;
        end
        
        function varargout = jointLimitConstraints(obj,varargin)
            varargout=cell(1,nargout);
            [varargout{:}] = jointLimitConstraints(obj.plant,varargin{:});
        end
        
        function index = getActuatedJoints(obj)
            index = getActuatedJoints(obj.plant);
        end
        
        function ptr = getMexModelPtr(obj)
            ptr = getMexModelPtr(obj.plant);
        end
        
        function [phi,Jphi] = closestDistance(obj,varargin)
            [phi,Jphi] = closestDistance(obj.plant,varargin{:});
        end
        
        function obj = addLinksToCollisionFilterGroup(obj,linknames,collision_fg_name,robotnums)
            obj.plant = addLinksToCollisionFilterGroup(obj.plant,linknames,collision_fg_name,robotnums);
        end
        
        function out = name(obj)
            out = obj.plant.name;
        end
        
        function fr = getParamFrame(model)
            fr = getParamFrame(model.plant);
        end
        
        function model = setParams(model,p)
            model.plant = setParams(model.plant,p);
        end
        
        function terrain_contact_point_struct = getTerrainContactPoints(obj,varargin)
            terrain_contact_point_struct = getTerrainContactPoints(obj.plant,varargin{:});
        end
        
        function varargout = terrainContactPositions(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = terrainContactPositions(obj.plant,varargin{:});
        end
        
        function varargout = terrainContactJacobianDotTimesV(obj,varargin)
            varargout = cell(1,nargout);
            [varargout{:}] = terrainContactJacobianDotTimesV(obj.plant,varargin{:});
        end
        
        function distance = collisionRaycast(obj, kinsol, origin, point_on_ray, use_margins)
            if nargin < 5
                use_margins = true;
            end
            distance = collisionRaycast(obj.plant, kinsol, origin, point_on_ray, use_margins);
        end
        
        
    end
    
end

