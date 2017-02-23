classdef VariationalTimeSteppingRigidBodyManipulator < DrakeSystem
  % A discrete time system which simulates the manipulator equations
  % with contact / limits resolved using a complementarity formulation
  % of contact and a variational integrator for the dynamics.
  
  % Note: Does not like joint friction being included in C. Edit URDF to
  % set friction constant to zero before using this class.

  properties (Access=protected)
    manip  % the CT manipulator
    sensor % additional TimeSteppingRigidBodySensors (beyond the sensors attached to manip)
    dirty=true;
    
    % cache stuff for warm starting
    cache_x
    cache_f
    cache_q1
    cache_q2
    cache_q3
  end

  properties (SetAccess=protected)
    timestep
    integrator
    twoD=false
    multiple_contacts = false;
    Nd = 2;
  end
  
  properties (Constant)
    MIDPOINT = 2;  % 2nd order variational midpoint integrator
    SIMPSON = 3;   % DEFAULT - 3rd order Simpson's rule variational integrator
  end

  methods
    function obj=VariationalTimeSteppingRigidBodyManipulator(manipulator_or_urdf_filename,timestep,options)
      if (nargin<3) options=struct(); end
      if ~isfield(options,'twoD') options.twoD = false; end
        
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

      if isfield(options, 'multiple_contacts')
        typecheck(options.multiple_contacts, 'logical');
        obj.multiple_contacts = options.multiple_contacts;
      end
      
      if isfield(options,'integration_method')
          obj.integrator = options.integration_method;
      else
          obj.integrator = VariationalTimeSteppingRigidBodyManipulator.SIMPSON;
      end

      obj.timestep = timestep;
      
      obj.cache_x = SharedDataHandle(0);
      obj.cache_f = SharedDataHandle(0);
      obj.cache_q1 = SharedDataHandle(0);
      obj.cache_q2 = SharedDataHandle(0);
      obj.cache_q3 = SharedDataHandle(0);

      obj = setSampleTime(obj,[timestep;0]);

      obj = compile(obj);
    end

    function checkDirty(obj)
      if (obj.dirty)
        error('You''ve changed something about this model and need to manually compile it.  Use obj=compile(obj).');
      end
    end

    function manip = getManipulator(obj)
      manip = obj.manip;
    end

    function y = output(obj,t,x,u)
      checkDirty(obj);

      if ~isDirectFeedthrough(obj)
        u=[];
      end
      if isa(obj.getStateFrame(),'MultiCoordinateFrame')
        x_manip = double(Point(obj.getStateFrame(),x).inFrame(obj.manip.getStateFrame()));
      else
        x_manip = x;
      end
      y = obj.manip.output(t,x_manip,u);
      for i=1:length(obj.sensor)
        y = [y; obj.sensor{i}.output(obj,i+1,t,x,u)];
      end
    end

    function model = compile(model)
      w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
      model.manip = model.manip.compile();
      warning(w);

      model = setNumDiscStates(model,model.manip.getNumContStates());
      model = setNumInputs(model,model.manip.getNumInputs());

%       if (model.position_control)
%         index = getActuatedJoints(model.manip);
%         pdff = pdcontrol(model,eye(model.num_u),eye(model.num_u));
%         model = setInputLimits(model,model.manip.joint_limit_min(index),model.manip.joint_limit_max(index));
%         model = setInputFrame(model,getInputFrame(pdff));
%       else
        model = setInputLimits(model,model.manip.umin,model.manip.umax);
        model = setInputFrame(model,getInputFrame(model.manip));
%       end
      model = setStateFrame(model,getStateFrame(model.manip));

      if ~isempty(model.sensor)
        feedthrough = model.manip.isDirectFeedthrough;
        outframe{1} = getOutputFrame(model.manip);
        stateframe{1} = getStateFrame(model.manip);
        for i=1:length(model.sensor)
          model.sensor{i} = model.sensor{i}.compile(model,model.manip);
          outframe{i+1} = model.sensor{i}.constructFrame(model);
          if isa(model.sensor{i},'TimeSteppingRigidBodySensorWithState')
            stateframe{i+1} = model.sensor{i}.constructStateFrame(model);
          end
          feedthrough = feedthrough || model.sensor{i}.isDirectFeedthrough;
        end
        fr = MultiCoordinateFrame.constructFrame(outframe);
        state_fr = MultiCoordinateFrame.constructFrame(stateframe);
        if ~isequal_modulo_transforms(fr,getOutputFrame(model))
          model = setNumOutputs(model,fr.dim);
          model = setOutputFrame(model,fr);
        end
        if ~isequal_modulo_transforms(state_fr,getStateFrame(model))
          model = setNumDiscStates(model,state_fr.dim);
          model = setStateFrame(model,state_fr);
        end
        model = setDirectFeedthrough(model,feedthrough);
      else
        model = setNumOutputs(model,getNumOutputs(model.manip));
        model = setOutputFrame(model,getOutputFrame(model.manip));
        model = setDirectFeedthrough(model,model.manip.isDirectFeedthrough);
      end
      model.dirty = false;
    end

    function x0 = getInitialState(obj)
      if ~isempty(obj.initial_state)
        x0 = obj.initial_state;
        return;
      end
      
      x0 = obj.manip.getInitialState();
      for i=1:length(obj.sensor)
        if isa(obj.sensor{i},'TimeSteppingRigidBodySensorWithState')
          x0 = [x0; obj.sensor{i}.getInitialState(obj)];
        end
      end
    end
    
    function [xdn, dxdn] = update(obj,t,x,u)
        
        h = obj.timestep;
        Nq = obj.manip.num_positions;
        Nd = obj.Nd;
        
        %vToqdot = obj.manip.vToqdot(kinsol); %not sure what this is for yet...
        zguess = zeros(2*Nq+2+2*Nd,1);
        if (length(obj.cache_x.data) == length(x)) && all(obj.cache_x.data == x)
            %we've been here before...
            q1 = x(1:Nq);
            v1 = obj.cache_x.data((Nq+1):end);
            p1 = obj.DLT(obj.cache_q1.data,obj.cache_q2.data,obj.cache_q3.data);
            zguess = [q1+(h/2)*v1; q1+h*v1; obj.cache_f.data];
        else
            %Set up intitial conditions
            q1 = x(1:Nq);
            v1 = x((Nq+1):end);
            M = manipulatorDynamics(obj.manip, q1, v1);
            p1 = M*v1;
            zguess = [q1+(h/2)*v1; q1+h*v1; zeros(2+2*Nd,1)];
        end
        
        %Solve NLP
        fmincon_opts = optimoptions('fmincon','Algorithm','sqp','Display','off','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'OptimalityTolerance',1e-8,'ConstraintTolerance',1e-8);
        zopt = fmincon(@(z)SimpsonObjFun(obj,q1,z), zguess, [], [], [], [], [], [], @(z)SimpsonConFun(obj,p1,q1,z), fmincon_opts);
        
        %update state
        q2 = zopt(1:Nq);
        q3 = zopt(Nq+(1:Nq));
        v3 = (q1 - 4*q2 + 3*q3)/h;
        xdn = [q3; v3];
        
        %update cache
        obj.cache_x.data = xdn;
        obj.cache_f.data = zopt((2*Nq+1):end);
        obj.cache_q1.data = q1;
        obj.cache_q2.data = q2;
        obj.cache_q3.data = q3;
        
        for i=1:length(obj.sensor)
            if isa(obj.sensor{i},'TimeSteppingRigidBodySensorWithState')
                if (nargout>1)
                    [obj,xdn_sensor,df_sensor] = update(obj.sensor{i},obj,t,x,u);
                else
                    [obj,xdn_sensor] = update(obj.sensor{i},obj,t,x,u);
                end
                xdn = [xdn;xdn_sensor];
                if (nargout>1)
                    df = [df; df_sensor];
                end
            end
        end
        
    end
    
    function [e, de] = SimpsonObjFun(obj, q1, z)
        %Minimize kinetic energy at end of time step
        Nq = length(q1);
        Nd = obj.Nd;
        Nx = length(z)-(2+2*Nd);
        
        h = obj.timestep;
        q2 = z(1:Nq);
        q3 = z(Nq+(1:Nq));
        
        v3 = (q1 - 4*q2 + 3*q3)/h;
        dv3 = [(-4/h)*eye(Nq), (3/h)*eye(Nq)];
        
        [M,~,~,dM] = manipulatorDynamics(obj.manip, q3, v3);
        dM = reshape(dM,Nq*Nq,Nx);
        dMdq = dM(:,1:Nq);
        C5 = reshape(dMdq*v3,Nq,Nq); %this is just the matrix C(q,qd)
        
        e = 0.5*v3'*M*v3;
        de = [v3'*M*dv3 + [zeros(1,Nq), 0.5*v3'*C5], zeros(1, 2+2*Nd)];
    end
    
    function [c, ceq, dc, dceq] = SimpsonConFun(obj, p1, q1, z)
        %This is repackaged for fmincon
        [c, dc] = SimpsonIneqCon(obj, q1, z);
        c = -c;
        dc = -dc';
        [ceq, dceq] = SimpsonEqCon(obj, p1, q1, z);
        dceq = dceq';
    end
    
    function [r, dr] = SimpsonEqCon(obj, p1, q1, z)
        Nq = length(q1);
        Nd = obj.Nd;
        h = obj.timestep;
        q2 = z(1:Nq);
        q3 = z(Nq+(1:Nq));
        
        %Contact force coefficients
        c1 = z(2*Nq+1);
        b1 = z(2*Nq+1+(1:Nd));
        c2 = z(2*Nq+1+Nd+1);
        b2 = z(2*Nq+1+Nd+1+(1:Nd));
        
        [r, dr] = SimpsonDEL(obj,p1,q1,q2,q3);
        
        %Add in contact forces
        kin1 = obj.manip.doKinematics(q1);
        [~,~,~,~,~,~,~,~,n1,D1] = obj.manip.contactConstraints(kin1, obj.multiple_contacts);
        n1 = n1';
        D1 = cell2mat(D1)';
        kin2 = obj.manip.doKinematics(q2);
        [~,~,~,~,~,~,~,~,n2,D2] = obj.manip.contactConstraints(kin2, obj.multiple_contacts);
        n2 = n2';
        D2 = cell2mat(D2)';
        
        r = r + [(h/3)*(n1*c1 + D1*b1); (2*h/3)*(n2*c2 + D2*b2)];
        dr = [dr, [(h/3)*n1, (h/3)*D1, zeros(Nq,1+Nd); zeros(Nq,1+Nd), (2*h/3)*n2, (2*h/3)*D2]];
    end
    
    function [p, dp] = SimpsonIneqCon(obj, q1, z)
        Nq = length(q1);
        Nd = obj.Nd;
        h = obj.timestep;
        q2 = z(1:Nq);
        q3 = z(Nq+(1:Nq));
        
        %Contact force coefficients
        c1 = z(2*Nq+1);
        b1 = z(2*Nq+1+(1:Nd));
        c2 = z(2*Nq+1+Nd+1);
        b2 = z(2*Nq+1+Nd+1+(1:Nd));
        
        %Contact force basis vectors
        kin2 = obj.manip.doKinematics(q2);
        [phi2,~,~,~,~,~,~,mu,n2,D2] = obj.manip.contactConstraints(kin2, obj.multiple_contacts);
        n2 = n2';
        D2 = cell2mat(D2)';
        kin3 = obj.manip.doKinematics(q3);
        [phi3,~,~,~,~,~,~,~,n3,D3] = obj.manip.contactConstraints(kin3, obj.multiple_contacts);
        n3 = n3';
        D3 = cell2mat(D3)';
        
        %velocity at last knot point
        v3 = (q1 - 4*q2 + 3*q3)/h;
        dv3 = [(-4/h)*eye(Nq), (3/h)*eye(Nq)];
        
        p = [phi2; phi3; c1; b1; c2; b2; %height + forces have to be positive
            -phi3*[c1; b1; c2; b2]; %contact forces can only if in contact at end of time step
            -(n3'*v3)*[c1; b1; c2; b2]; %contact forces can only act if normal velocity component <= 0 at end of time step
             mu*c1-b1; mu*c2-b2]; %friction cone
          
        dp = [n2', zeros(1,Nq+2+2*Nd); %height has to be positive
              zeros(1,Nq), n3', zeros(1, 2+2*Nd); %height has to be positive
              zeros(2+2*Nd, 2*Nq), eye(2+2*Nd); %forces have to be positive
              zeros(2+2*Nd,Nq), -[c1; b1; c2; b2]*n3', -phi3*eye(2+2*Nd); %contact forces can only if in contact at end of time step
              -[c1; b1; c2; b2]*(n3'*dv3), -(n3'*v3)*eye(2+2*Nd); %contact forces can only act if normal velocity component <= 0 at end of time step
              zeros(Nd,2*Nq), mu*ones(Nd,1), -eye(Nd), zeros(Nd,1+Nd); %friction cone
              zeros(Nd,2*Nq), zeros(Nd, 1+Nd), mu*ones(Nd,1), -eye(Nd)]; %friction cone
    end
    
    function [r, dr] = MidpointDEL(obj, q1, q2, q3)
        r = (h/2)*obj.D1L((q1+q2)/2,(q2-q1)/h) + obj.D2L((q1+q2)/2,(q2-q1)/h) + (h/2)*obj.D1L((q2+q3)/2,(q3-q2)/h) - obj.D2L((q2+q3)/2,(q3-q2)/h);
        dr = (h/4)*obj.D1D1L((q2+q3)/2,(q3-q2)/h) + (1/2)*obj.D1D2L((q2+q3)/2,(q3-q2)/h) - (1/2)*obj.D1D2L((q2+q3)/2,(q3-q2)/h) - (1/h)*obj.D2D2L((q2+q3)/2,(q3-q2)/h);
    end
    
    function [r, dr] = SimpsonDEL(obj, p1, q1, q2, q3)
        h = obj.timestep;
        
        %velocities at timestep k+1 knot points
        v1 = (-3*q1 + 4*q2 - q3)/h;
        v2 = (q3 - q1)/h;
        v3 = (q1 - 4*q2 + 3*q3)/h;
        
        r = [p1 + (h/6)*obj.D1L(q1,v1) - (1/2)*obj.D2L(q1,v1) - (2/3)*obj.D2L(q2,v2) + (1/6)*obj.D2L(q3,v3);
             (2/3)*obj.D2L(q1,v1) + (2*h/3)*obj.D1L(q2,v2) - (2/3)*obj.D2L(q3,v3)];
        
        dr = [(2/3)*obj.D1D2L(q1,v1) - (2/h)*obj.D2D2L(q1,v1) - (2/3)*obj.D1D2L(q2,v2) - (2/(3*h))*obj.D2D2L(q3,v3), -(1/6)*obj.D1D2L(q1,v1) + (1/(2*h))*obj.D2D2L(q1,v1) - (2/(3*h))*obj.D2D2L(q2,v2) + (1/6)*obj.D1D2L(q3,v3) + (1/(2*h))*obj.D2D2L(q3,v3);
              (8/(3*h))*obj.D2D2L(q1,v1) + (2*h/3)*obj.D1D1L(q2,v2) + (8/(3*h))*obj.D2D2L(q3,v3), -(2/(3*h))*obj.D2D2L(q1,v1) + (2/3)*obj.D1D2L(q2,v2) - (2/3)*obj.D1D2L(q3,v3) - (2/h)*obj.D2D2L(q3,v3)];
    end
    
    function p = DLT(obj, q1, q2, q3)
        h = obj.timestep;
        
        v1 = (-3*q1 + 4*q2 - q3)/h;
        v2 = (q3 - q1)/h;
        v3 = (q1 - 4*q2 + 3*q3)/h;
        
        p = -(1/6)*obj.D2L(q1,v1) + (2/3)*obj.D2L(q2,v2) + (h/6)*obj.D1L(q3,v3) + (1/2)*obj.D2L(q3,v3);
    end
    
    function dL = D1L(obj,q,v)
        [~,Cg] = manipulatorDynamics(obj.manip, q, v);
        [~,G] = manipulatorDynamics(obj.manip, q, zeros(size(v)));
        
        Cv = Cg - G;
        dL = 0.5*Cv - G;
    end
    
    function dL = D2L(obj,q,v)
         M = manipulatorDynamics(obj.manip, q, v);
         dL = M*v;
    end
    
    function d2L = D1D1L(obj,q,v)
        Nq = length(q);
        Nv = length(v);
        
        [~,~,~,~,dCg] = manipulatorDynamics(obj.manip, q, v);
        [~,~,~,~,dG] = manipulatorDynamics(obj.manip, q, zeros(Nv,1));
        
        dCgdq = dCg(:,1:Nq);
        dGdq = dG(:,1:Nq);
        
        d2L = 0.5*(dCgdq - dGdq) - dGdq;
    end
    
    function d2L = D1D2L(obj,q,v)
        Nq = length(q);
        Nv = length(v);
        
        [~,~,~,dM] = manipulatorDynamics(obj.manip, q, v);
        dM = reshape(dM,Nq*Nq,Nq+Nv);
        dMdq = dM(:,1:Nq);
        d2L = reshape(dMdq*v,Nq,Nq); %this is just the matrix C(q,qd)
        
        %This should give the same answer...
        %Nq = length(q);
        %Nv = length(v);
        %[~,~,~,~,dCg] = manipulatorDynamics(obj.manip, q, v);
        %[~,~,~,~,dG] = manipulatorDynamics(obj.manip, q, zeros(Nv,1));
        %dCv = dCg - dG;
        %C = dCv(:,Nq+(1:Nv));
    end
    
    function d2L = D2D2L(obj,q,v)
        M = manipulatorDynamics(obj.manip, q, v);
        d2L = M;
    end
    
    
    function obj = addSensor(obj,sensor)
      if isa(sensor,'RigidBodySensor')
        obj.manip = obj.manip.addSensor(sensor);
      else
        typecheck(sensor,'TimeSteppingRigidBodySensor');
        obj.sensor{end+1} = sensor;
      end
    end

    function [obj,frame_id] = addFrame(obj,frame)
      [obj.manip,frame_id] = obj.manip.addFrame(frame);
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
