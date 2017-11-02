classdef LittleDog < RigidBodyManipulator %& LeggedRobot
  
  methods
    function obj = LittleDog(options)
      options.floating = true;
      if ~isfield(options,'terrain')
       options.terrain = RigidBodyFlatTerrain();          
      end
      w = warning('off','Drake:RigidBodyManipulator:UnsupportedContactPoints');
      obj = obj@RigidBodyManipulator('LittleDog.urdf',options);
      warning(w);
      
%      obj = addFoot(obj, 'front_left', 'front_left_foot_center');
%      obj = addFoot(obj, 'front_right', 'front_right_foot_center');
%      obj = addFoot(obj, 'back_left', 'back_left_foot_center');
%      obj = addFoot(obj, 'back_right', 'back_right_foot_center');
    end
    
    function x0 = home(obj)
      hip_roll = .1;
      hip_pitch = 1;
      knee = 1.55;
      x0 = Point(getStateFrame(obj));
      x0.front_right_hip_roll = -hip_roll;
      x0.front_right_hip_pitch = hip_pitch;
      x0.front_right_knee = -knee;
      x0.front_left_hip_roll = hip_roll;
      x0.front_left_hip_pitch = hip_pitch;
      x0.front_left_knee = -knee;
      x0.back_right_hip_roll = -hip_roll;
      x0.back_right_hip_pitch = -hip_pitch;
      x0.back_right_knee = knee;
      x0.back_left_hip_roll = hip_roll;
      x0.back_left_hip_pitch = -hip_pitch;
      x0.back_left_knee = knee;
%      x0 = resolveConstraints(obj,x0);  % andres says calling this results
%      in snopt be slow thereafter.  In fact, it only solves for one var:
      x0.base_z = 0.146;
    end
    
    function xstar = defaultFixedPoint(obj)
      xstar = home(obj); %findFixedPoint(obj,home(obj));
    end
    
    function plan = planFootsteps(obj)
      error('not implemented yet')
    end
    
  end
  
  methods (Static=true)
    function runPassive
      r = LittleDog;
      sys = TimeSteppingRigidBodyManipulator(r,0.01);
      x0 = home(r);
      
      v = constructVisualizer(sys);
      xtraj = simulate(sys,[0 2],double(x0));
      playback(v,xtraj);
    end
    
    function runPDhome
      r = LittleDog;
      sys = TimeSteppingRigidBodyManipulator(r,0.001);
      v = constructVisualizer(sys);
      
      % construct PD control 
      Kp = 2*eye(12);
      Kd = diag([0.5; 0.5; 0.16; 0.5; 0.5; 0.16; 0.5; 0.5; 0.16; 0.5; 0.5; 0.16]);
      sys = pdcontrol(sys,Kp,Kd);
      
      x0 = home(r);
      qa0 = x0(getActuatedJoints(r));

      % send in constant position reference
      sys = cascade(setOutputFrame(ConstOrPassthroughSystem(qa0),getInputFrame(sys)),sys);

      if (0)
        sys = cascade(sys,v);
        simulate(sys,[0 3],double(x0));
      else
        xtraj = simulate(sys,[0 3],double(x0));
        playback(v,xtraj);
      end
    end
    
    function [e, rmse, esamp] = runPDTraj(xtraj)
      r = LittleDog;
      sys = TimeSteppingRigidBodyManipulator(r,0.0002);
      v = constructVisualizer(sys);
      
      % construct PD control 
      Kp = 3*eye(12);
      Kd = .2*diag([1; 1; .5; 1; 1; .5; 1; 1; .5; 1; 1; .5]);
      
      sys = pdcontrol(sys,Kp,Kd);
      
      qatraj = xtraj(getActuatedJoints(r));

      % send in joint angle reference
      sys = cascade(setOutputFrame(qatraj,getInputFrame(sys)),sys);
      
      x0 = xtraj.eval(0);
      phi = r.contactConstraints(x0(1:18));
      while min(phi) < 0
        x0(3) = x0(3) + 0.01;
        phi = r.contactConstraints(x0(1:18));
      end

      cltraj = simulate(sys,[xtraj.tspan(1), xtraj.tspan(2)+.05],x0);
      playback(v,cltraj,struct('slider',true));
      
      tsamp = floor(1000*linspace(cltraj.tspan(1),cltraj.tspan(2),100))/1000;
      esamp = zeros(36,100);
      rmse = zeros(100,1);
      for k = 1:length(tsamp)
          if tsamp(k) <= xtraj.tspan(2)
              xk = xtraj.eval(tsamp(k));
          else
              xk = xtraj.eval(xtraj.tspan(2));
          end
          esamp(:,k) = xk - cltraj.eval(tsamp(k));
          rmse(k) = sqrt(esamp(:,k)'*esamp(:,k))/36;
      end
      e = sum(rmse)/100;
      
    end
    
    function runOLTraj(xtraj,utraj)
      r = LittleDog;
      sys = TimeSteppingRigidBodyManipulator(r,0.0002);
      v = constructVisualizer(sys);

      sys = cascade(setOutputFrame(utraj,getInputFrame(sys)),sys);
      
      x0 = xtraj.eval(0);
      phi = r.contactConstraints(x0(1:18));
      while min(phi) < 0
        x0(3) = x0(3) + 0.01;
        phi = r.contactConstraints(x0(1:18));
      end

      cltraj = simulate(sys,[xtraj.tspan(1), xtraj.tspan(2)+.1],x0);
      playback(v,cltraj,struct('slider',true));
    end
  end
end
