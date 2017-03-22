classdef VariationalTimeSteppingRigidBodyManipulator < DrakeSystem
  % A discrete time system which simulates the manipulator equations
  % with contact / limits resolved using a complementarity formulation
  % of contact and a variational integrator for the dynamics.
  
  % Note: Does not like joint friction being included in C. Edit URDF to
  % set friction constant to zero before using this class.

  properties %(Access=protected)
    manip  % the CT manipulator
    sensor % additional TimeSteppingRigidBodySensors (beyond the sensors attached to manip)
    dirty=true;
    
    % cache stuff for warm starting
    cache_x
    cache_q1
    cache_q2
    cache_q3
  end

  properties (SetAccess=protected)
    timestep
    tolerance = 1e-7;
    damping = 1e-3;
    sdamping = .1;
    integrator
    twoD=false
    multiple_contacts = false;
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
        Nq = obj.getNumPositions;
        
        if (length(obj.cache_x.data) == length(x)) && all(obj.cache_x.data == x)
            %we've been here before...
            q1 = x(1:Nq);
            v1 = obj.cache_x.data((Nq+1):end);
            p1 = obj.DLT(obj.cache_q1.data,obj.cache_q2.data,obj.cache_q3.data);
        else
            %Set up intitial conditions
            q1 = x(1:Nq);
            v1 = x((Nq+1):end);
            M = manipulatorDynamics(obj.manip, q1, v1);
            p1 = M*v1;
        end
        
        kin1 = obj.manip.doKinematics(q1);
        colopts = struct();
        [phi1,~,~,~,~,~,~,~,n1,D1] = obj.manip.contactConstraints(kin1, obj.multiple_contacts, colopts);
        
        Np = size(n1,1);
        Nd = length(D1);
        
        D1 = reshape(cell2mat(D1')',Nq,Np*Nd)';
        
        zguess = [q1; q1; zeros(2*Np+2*Np*Nd+2*Np,1); 1];
        
        %Solve Complementarity Problem
        if Np == 0
            %No contact -- yay!
            q2 = zguess(1:Nq);
            q3 = zguess(Nq+(1:Nq));
            r = 1;
            while max(abs(r)) > obj.tolerance
                [r, dr] = obj.SimpsonDEL(p1, q1, q2, q3);
                delta = -dr\r;
                q2 = q2 + delta(1:Nq);
                q3 = q3 + delta(Nq+(1:Nq));
            end
        else
            zopt = obj.dnsolver(q1,p1,zguess,Np,Nd);
            q2 = zopt(1:Nq);
            q3 = zopt(Nq+(1:Nq));
        end
        
        %update state
        v3 = (q1 - 4*q2 + 3*q3)/h;
        xdn = [q3; v3];
        
        %update cache
        obj.cache_x.data = xdn;
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
    
    function z = dnsolver(obj,q1,p1,zguess,Np,Nd)
        tol = obj.tolerance;
        Nq = length(q1);
        
        z = zguess;
        f = 1;
        L = diag([zeros(2*Nq,1); obj.damping*ones(2*(Np+Np*Nd),1); zeros(2*Np,1); obj.sdamping]);
        while max(abs(f)) > tol
            
            f = obj.SimpsonResidual(q1,p1,z,Np,Nd);
            J = f'*f;
            
            %--------- First Derivatives ---------%
            dz = .5*tol*eye(length(z));
            df = zeros(length(f), length(z));
            for k = 1:length(z)
                df(:,k) = (obj.SimpsonResidual(q1,p1,z+dz(:,k),Np,Nd)-obj.SimpsonResidual(q1,p1,z-dz(:,k),Np,Nd))/tol;
            end
            %-------------------------------------%
            
            %Damped + Regularized Newton descent direction
            [Q,R] = qr([df; L],0);
            delta = -R\(Q(1:length(f),:)'*f);
            
            znew = z+delta;
            
            fnew = obj.SimpsonResidual(q1,p1,znew,Np,Nd);
            Jnew = fnew'*fnew;
            
            if Jnew > J
                alpha = (1/2);
                while Jnew > J && alpha > tol
                    znew = z + alpha*delta;
                    fnew = obj.SimpsonResidual(q1,p1,znew,Np,Nd);
                    Jnew = fnew'*fnew;
                    alpha = alpha/2;
                end
                if alpha < tol
                    znew = z;
                end
            end
            z = znew;
        end
    end
    
    function f = SimpsonResidual(obj,q1,p1,z,Np,Nd)
        mu = 1; %This is currently hard coded in Drake...
        
        Nq = length(q1);

        h = obj.timestep;
        
        %z vector is stacked [q_{1/2}; q_2; c1; b1; c2; b2; eta; psi; s]
        
        %Configurations
        q2 = z(1:Nq);
        q3 = z(Nq+(1:Nq));
        
        %Contact force coefficients
        c1 = z(2*Nq+(1:Np));
        b1 = z(2*Nq+Np+(1:Np*Nd));
        c2 = z(2*Nq+Np+Np*Nd+(1:Np));
        b2 = z(2*Nq+2*Np+Np*Nd+(1:Np*Nd));
        
        %Normal contact velocity
        eta = z(2*Nq+2*Np+2*Np*Nd+(1:Np));
        
        %Tangential contact velocity
        psi = z(2*Nq+2*Np+2*Np*Nd+Np+(1:Np));
        
        s = z(end);
        
        %velocity at last knot point
        v3 = (q1 - 4*q2 + 3*q3)/h;
%         dv3 = [(-4/h)*eye(Nq), (3/h)*eye(Nq)];
        
        %Get contact basis
        e = ones(Nd,1);
        
        kin1 = obj.manip.doKinematics(q1);
        [phi1,~,~,~,~,~,~,~,n1,D1] = obj.manip.contactConstraints(kin1, obj.multiple_contacts);
        
        kin2 = obj.manip.doKinematics(q2);
        [phi2,~,~,~,~,~,~,~,n2,D2] = obj.manip.contactConstraints(kin2, obj.multiple_contacts);
        
        kin3 = obj.manip.doKinematics(q3);
        [phi3,~,~,~,~,~,~,~,n3,D3] = obj.manip.contactConstraints(kin3, obj.multiple_contacts);
        
        D1 = reshape(cell2mat(D1')',Nq,Np*Nd)';
        D2 = reshape(cell2mat(D2')',Nq,Np*Nd)';
        D3 = reshape(cell2mat(D3')',Nq,Np*Nd)';
            
        r_del = obj.SimpsonDEL(p1, q1, q2, q3);
        
        r_f = r_del + [(h/3)*(n1'*c1 + D1'*b1); (2*h/3)*(n2'*c2 + D2'*b2)];
%         dr_f = [dr_del, [(h/3)*n1', (h/3)*D1', zeros(Nq,Np+Np*Nd); zeros(Nq,Np+Np*Nd), (2*h/3)*n2', (2*h/3)*D2'], zeros(2*Nq,2*Np+1)];
        
        [fb1, dfba1, dfbb1, dfbs1] = obj.smoothFB(phi3, c1, s);
        [fb2, dfba2, dfbb2, dfbs2] = obj.smoothFB(kron(phi3,ones(Nd,1)), b1, s);
        [fb3, dfba3, dfbb3, dfbs3] = obj.smoothFB(phi3, c2, s);
        [fb4, dfba4, dfbb4, dfbs4] = obj.smoothFB(kron(phi3,ones(Nd,1)), b2, s);
        
        [fb5, dfba5, dfbb5, dfbs5] = obj.smoothFB(eta, c1, s);
        [fb6, dfba6, dfbb6, dfbs6] = obj.smoothFB(kron(eta,ones(Nd,1)), b1, s);
        [fb7, dfba7, dfbb7, dfbs7] = obj.smoothFB(eta, c2, s);
        [fb8, dfba8, dfbb8, dfbs8] = obj.smoothFB(kron(eta,ones(Nd,1)), b2, s);
        
        [fb9, dfba9, dfbb9, dfbs9] = obj.smoothFB(eta, (eta - n3*v3), s);
        [fb10, dfba10, dfbb10, dfbs10] = obj.smoothFB((mu*c1 - kron(eye(Np),e')*b1), psi, s);
        [fb11, dfba11, dfbb11, dfbs11] = obj.smoothFB((mu*c2 - kron(eye(Np),e')*b2), psi, s);
        [fb12, dfba12, dfbb12, dfbs12] = obj.smoothFB((kron(psi,e) + D3*v3), b1, s);
        [fb13, dfba13, dfbb13, dfbs13] = obj.smoothFB((kron(psi,e) + D3*v3), b2, s);
        
        f = [r_f; fb1; fb2; fb3; fb4; fb5; fb6; fb7; fb8; fb9; fb10; fb11; fb12; fb13; exp(s)-1];
%         df = [dr_f;
%               zeros(1,2*Nq+2*Np*Nd+4*Np), exp(s);
%               zeros(2*(Np+Np*Nd),Nq), dfba1*kron(n3,ones(2*(Nd+1),1)), dfbb1, zeros(2*(Np+Np*Nd),2*Np), dfbs1;
%               zeros(2*(Np+Np*Nd),2*Nq), dfbb2, dfba2*kron(eye(Np),ones(2*(Nd+1),1)), zeros(2*(Np+Np*Nd),Np), dfbs2;
%               -dfbb3*n3*dv3, zeros(Np,2*(Np+Np*Nd)), dfba3+dfbb3, zeros(Np), dfbs3;
%               zeros(Np,2*Nq), mu*dfba4, -dfba4*kron(eye(Np),e'), zeros(Np,Np+Np*Nd), zeros(Np), dfbb4, dfbs4;
%               zeros(Np,2*Nq), zeros(Np,Np+Np*Nd), mu*dfba5, -dfba5*kron(eye(Np),e'), zeros(Np), dfbb5, dfbs5;
%               dfba6*D3*dv3, zeros(Np*Nd,Np), dfbb6, zeros(Np*Nd,Np+Np*Nd), zeros(Np*Nd,Np), dfba6*kron(eye(Np),e), dfbs6;
%               dfba7*D3*dv3, zeros(Np*Nd,Np+Np*Nd), zeros(Np*Nd,Np), dfbb7, zeros(Np*Nd,Np), dfba7*kron(eye(Np),e), dfbs7];
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
    
    function [r, dr] = MidpointDEL(obj, q1, q2, q3)
        r = (h/2)*obj.D1L((q1+q2)/2,(q2-q1)/h) + obj.D2L((q1+q2)/2,(q2-q1)/h) + (h/2)*obj.D1L((q2+q3)/2,(q3-q2)/h) - obj.D2L((q2+q3)/2,(q3-q2)/h);
        dr = (h/4)*obj.D1D1L((q2+q3)/2,(q3-q2)/h) + (1/2)*obj.D1D2L((q2+q3)/2,(q3-q2)/h) - (1/2)*obj.D1D2L((q2+q3)/2,(q3-q2)/h) - (1/h)*obj.D2D2L((q2+q3)/2,(q3-q2)/h);
    end
    
    function r = SimpsonDEL(obj, p1, q1, q2, q3)
        h = obj.timestep;
        
        %velocities at timestep k+1 knot points
        v1 = (-3*q1 + 4*q2 - q3)/h;
        v2 = (q3 - q1)/h;
        v3 = (q1 - 4*q2 + 3*q3)/h;
        
        [d1L_1, d2L_1] = obj.SimpsonLDerivs(q1,v1);
        [d1L_2, d2L_2] = obj.SimpsonLDerivs(q2,v2);
        [d1L_3, d2L_3] = obj.SimpsonLDerivs(q3,v3);
        
%         [d1L_1, d2L_1, d1d1L_1, d1d2L_1, d2d2L_1] = obj.SimpsonLDerivs(q1,v1);
%         [d1L_2, d2L_2, d1d1L_2, d1d2L_2, d2d2L_2] = obj.SimpsonLDerivs(q2,v2);
%         [d1L_3, d2L_3, d1d1L_3, d1d2L_3, d2d2L_3] = obj.SimpsonLDerivs(q3,v3);

        r = [p1 + (h/6)*d1L_1 - (1/2)*d2L_1 - (2/3)*d2L_2 + (1/6)*d2L_3;
             (2/3)*d2L_1 + (2*h/3)*d1L_2 - (2/3)*d2L_3];
        
%         dr = [(2/3)*d1d2L_1 - (2/h)*d2d2L_1 - (2/3)*d1d2L_2 - (2/(3*h))*d2d2L_3, -(1/6)*d1d2L_1 + (1/(2*h))*d2d2L_1 - (2/(3*h))*d2d2L_2 + (1/6)*d1d2L_3 + (1/(2*h))*d2d2L_3;
%               (8/(3*h))*d2d2L_1 + (2*h/3)*d1d1L_2 + (8/(3*h))*d2d2L_3, -(2/(3*h))*d2d2L_1 + (2/3)*d1d2L_2 - (2/3)*d1d2L_3 - (2/h)*d2d2L_3];
    end
    
    function p = DLT(obj, q1, q2, q3)
        h = obj.timestep;
        
        v1 = (-3*q1 + 4*q2 - q3)/h;
        v2 = (q3 - q1)/h;
        v3 = (q1 - 4*q2 + 3*q3)/h;
        
        p = -(1/6)*obj.D2L(q1,v1) + (2/3)*obj.D2L(q2,v2) + (h/6)*obj.D1L(q3,v3) + (1/2)*obj.D2L(q3,v3);
    end
    
    function dL = D1L(obj,q,v)
        Nq = length(q);
        Nv = length(v);
        [~,G,~,dM] = manipulatorDynamics(obj.manip, q, zeros(Nv,1));
        dM = reshape(dM,Nq*Nq,Nq+Nv);
        dMdq = dM(:,1:Nq);
        dL = 0.5*dMdq'*kron(v,v) - G;
    end
    
    function dL = D2L(obj,q,v)
        M = manipulatorDynamics(obj.manip, q, v);
        dL = M*v;
    end
    
    function [d1L, d2L] = SimpsonLDerivs(obj,q,v)
        Nq = length(q);
        Nv = length(v);
        
        [M,G,~,dMdqx,dG] = manipulatorDynamics(obj.manip, q, zeros(Nv,1));
        
        dMdqx = reshape(dMdqx,Nq*Nq,Nq+Nv);
        dM = dMdqx(:,1:Nq);
        
        d1L = 0.5*dM'*kron(v,v) - G;
        d2L = M*v;
        
%         d2M = zeros(Nq*Nq*Nq,Nq);
%         dq = 5e-8*eye(Nq);
%         for k = 1:Nq
%             [~,~,~,dMp] = manipulatorDynamics(obj.manip, q+dq(:,k), zeros(Nv,1));
%             [~,~,~,dMm] = manipulatorDynamics(obj.manip, q-dq(:,k), zeros(Nv,1));
%             d2M(:,k) = vec(dMp(:,1:Nq)-dMm(:,1:Nq))/1e-7;
%         end
%         
%         d1d1L = 0.5*kron(eye(Nq),kron(v',v'))*d2M - dG(:,1:Nq);
%         d1d2L = kron(v',eye(Nq))*dM;
%         d2d2L = M;
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
