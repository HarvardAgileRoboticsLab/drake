classdef IndirectTrajectoryOptimization
  %INDIRECTTRAJECTORYOPTIMIZATION An abstract class for indirect approaches
  % to trajectory optimization.
  %
  % Generally considers cost functions of the form:
  % int(f(x(t),u(t)) + g(T,xf)
  %
  % Subclasses must implement the abstract method:
  %  [xtraj,utraj,J,info] = solveTraj(obj,t_init,traj_init)
  %
  % This class assumes that there are a fixed number (N) time steps, and
  % that the trajectory is discreteized into timesteps h (N-1), state x
  % (N), and control input u (N)
  %
  % This implementation assumes that all constraints and costs are
  % time-invariant.

  properties (SetAccess = protected)
    N       % number of timesteps
    duration % the total duration
    options % options, yup
    plant   % the plant
    running_cost_function % function handle that evaluates running cost
    final_cost_function % function handle that evaluates final cost
  end

  methods
    function obj = IndirectTrajectoryOptimization(plant,N,duration,options)
    % function obj =
    % DirectTrajectoryOptimization(plant,initial_cost,running_cost,final_cost,...
    % t_init,traj_init,T_span,constraints, options)
    % Trajectory optimization constructor
    % @param plant
    % @param N the number of time samples
    % @param duration  The total time for the trajectory
    % @param options (optional)

      if nargin < 4
        options = struct();
      end

      if ~plant.isTI
        error('Drake:IndirectTrajectoryOptimization:UnsupportedPlant','Only time-invariant plants are currently supported');
      end

      obj.N = N;
      obj.duration = duration;
      obj.options = options;
      obj.plant = plant;

    end
    
    function N = getN(obj)
      N = obj.N;
    end

    function obj = addTrajectoryDisplayFunction(obj,display_fun)
      % add a dispay function that gets called on every iteration of the
      % algorithm
      % @param display_fun a function handle of the form displayFun(t,x,u)
      %       where t is a 1-by-N, x is n-by-N and u is m-by-N
      
      obj = addDisplayFunction(obj,@(z) display_fun(z(obj.h_inds),z(obj.x_inds),z(obj.u_inds)));
    end
    
    function obj = addRunningCost(obj,running_cost_function);
        % Adds an integrated cost to all time steps
        % this cost is assumed to be time-invariant
        % @param running_cost_function a function handle g(dt,x,u)
        
        obj.running_cost_function = running_cost_function;
    end

    function obj = addFinalCost(obj,final_cost_function)
      % adds a cost to the final state and total time
      % @param final_cost_function a function handle f(T,xf)

      obj.final_cost_function = final_cost_function;
    end

    function utraj = reconstructInputTrajectory(obj,u)
      % default behavior is to use first order holds, but this can be
      % re-implemented by a subclass.
      try
          t = linspace(0,obj.duration,obj.N);
      catch
          t = linspace(0,obj.duration(2)-obj.duration(1),obj.N);
      end
      t = t(1:end-1);
      if isempty(u)
        utraj=[];
      else
        utraj = PPTrajectory(foh(t,u));
        utraj = utraj.setOutputFrame(obj.plant.getInputFrame);
      end
    end

    function xtraj = reconstructStateTrajectory(obj,x)
      % default behavior is to use first order holds, but this can be
      % re-implemented by a subclass.
      try
          t = linspace(0,obj.duration,obj.N);
      catch
          t = linspace(0,obj.duration(2)-obj.duration(1),obj.N);
      end
      xtraj = PPTrajectory(foh(t,x));
      xtraj = xtraj.setOutputFrame(obj.plant.getStateFrame);
    end

    function u0 = extractFirstInput(obj)
      % When using trajectory optimization a part of a model-predictive
      % control system, we often only need to extract u(0).  This method
      % does exactly that.
      u0 = obj.u(:,1);
    end
  end
  
  methods (Access = protected)
      function [l,q,Q] = running_cost(obj,h,x,u)
          if nargout > 1
              [l,df,d2f] = geval(obj.running_cost_function,h,x,u);
              q = df(2:end); %TODO: minimum time stuff
              Q = d2f(2:end,2:end);
          else
              l = obj.running_cost_function(h,x,u);
          end
      end
      
      function [l,qx,Qxx] = final_cost(obj,T,x)
          if nargout > 1
              [l,df,d2f] = geval(obj.final_cost_function,T,x);
              qx = df(2:end); %TODO: minimum time stuff
              Qxx = d2f(2:end,2:end);
          else
              l = obj.final_cost_function(T,x);
          end
      end
  end

  methods(Abstract)
    % Solve the optimization problem and return resulting trajectory
    % @param traj_init (optional) a structure containing Trajectory
    % objects specifying the initial guess for the system inputs
    %    traj_init.u , traj_init.x, ...
    % @default small random numbers
    [xtraj,utraj,J,info] = solveTraj(obj,tf,traj_init)
  end
end