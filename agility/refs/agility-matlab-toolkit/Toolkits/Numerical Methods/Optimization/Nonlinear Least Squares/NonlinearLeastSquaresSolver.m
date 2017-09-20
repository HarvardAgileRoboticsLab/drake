%NONLINEARLEASTSQUARESSOLVER Nonlinear least squares solver object.
%
% Syntax:
%   solver = NonlinearLeastSquaresSolver(problem);
%
% Description:
%   The Levenberg-Marquardt algorithm is a damped least squares method used
%   to solve nonlinear least squares problems.
%
% Notes:
%   Use small value for tau (1e-6) if x0 is believed to be close
%   to x*, otherwise use something larger (1e-3 < tau < 1).

% TODO: Add weighting, important for inverse kinematics
% TODO: Add getActiveBounds method

% Copyright 2016-2017 Mikhail S. Jones

classdef NonlinearLeastSquaresSolver < handle %#codegen

  % PROTECTED PROPERTIES ==================================================
  properties %(SetAccess = protected)
    % Nonlinear least squares problem
    problem
    % Decision variable vector
    x
    % Levenberg-Marquardt parameter
    levenbergMarquardtParameter
    % Solver iteration limit stopping criteria
    iterationLimit
    % Termination tolerance on x
    xTolerance
    % Termination tolerance on the gradient norm
    gTolerance
    % Solver status
    status
    % Current iteration count
    iteration
  end % properties

  % PUBLIC METHODS ========================================================
  methods (Access = public)
    function obj = NonlinearLeastSquaresSolver(problem)
    %NONLINEARLEASTSQUARESSOLVER Construct nonlinear least sqaures solver object.

      % Validate problem attributes
      validateattributes(problem, ...
        {'NonlinearLeastSquaresProblem'}, ...
        {'scalar'}, ...
        'NonlinearLeastSquaresSolver', ...
        'problem');
      
      % Set problem
      obj.problem = problem;

      % Set solution vector
      obj.x = zeros(problem.getNumberOfVariables, 1);
      
      % Initialize solver parameters
      obj.levenbergMarquardtParameter = 1e-3;
      obj.iterationLimit = 1e3;
      obj.xTolerance = sqrt(eps);
      obj.gTolerance = sqrt(eps);
      obj.status = OptimizationStatus.ITERATE;
    end % NonlinearLeastSquaresSolver

    function lsqnonlin(obj)
    %LSQNONLIN Use built in MATLAB solver.
    
      % Unpack problem structure into local variables
      [r, lb, ub] = obj.problem.unpack;
      
      % Solve problem
      obj.x = lsqnonlin(@(x) r.evaluate(x), obj.x, lb, ub);
      
      % Parse status
      % TODO
    end % lsqnonlin
    
    function solve(obj)
    %SOLVE Solve a nonlinear least squares problem.
    %
    % Syntax:
    %   obj.solve;
    %
    % Description:
    %   Solves the nonlinear least squares problem using the Levenberg
    %   Marquardt algorithm modified to allow bounds by using gradient
    %   projection.

      % Initialize solver status
      obj.status = OptimizationStatus.ITERATE;

      % Initialize iteration counter
      obj.iteration = 0;

      % Unpack problem structure into local variables
      [r, lb, ub] = obj.problem.unpack;

      % Unpack problem size into local variables
      nv = obj.problem.getNumberOfVariables;

      % Project initial point onto bounds
      [obj.x, r_x, J_x] = obj.projectOntoBounds(...
        r, obj.x, lb + 1e-9, ub - 1e-9);

      % Construct identity matrix used for computing step direction
      I = eye(nv);

      % Compute scalar objective to minimize (least squares)
      f_x = 0.5*(r_x'*r_x);

      % Compute gradient
      g_x = J_x'*r_x;

      % Check if norm of gradient is below tolerance
      if norm(g_x, Inf) <= obj.gTolerance
        % Set status and exit function
        obj.status = OptimizationStatus.G_TOLERANCE; return;
      end % if

      % Check if norm of gradient is below tolerance
      if isnan(norm(g_x, Inf))
        % Set status and exit function
        obj.status = OptimizationStatus.SINGULAR; return;
      end % if
      
      % Compute Gauss-Newton approximation of Hessian
      H_x = J_x'*J_x;

      % Set Levenberg-Marquardt parameters
      lambda = obj.levenbergMarquardtParameter*max(diag(H_x));
      nu = 2;

      % Main iteration loop
      while obj.status == OptimizationStatus.ITERATE
        % Advance iteration counter
        obj.iteration = obj.iteration + 1;
        
        % Check if iteration limit is exceeded
        if obj.iteration > obj.iterationLimit
          % Set status
          obj.status = OptimizationStatus.ITERATION_LIMIT;
        end % if

        % Compute step direction
        d = (H_x + lambda*I)\(-g_x);

        % Check if change of iterate is below tolerance
        if norm(d) <= obj.xTolerance*(obj.xTolerance + norm(obj.x))
          % Set status and exit loop
          obj.status = OptimizationStatus.X_TOLERANCE; break;
        end % if

        % Compute projected residual and Jacobian
        [x_new, r_new, J_new] = obj.projectOntoBounds(r, obj.x + d, lb, ub);

        % Compute scalar objective to minimize (least squares)
        f_new = 0.5*(r_new'*r_new);

        % Compute gain ratio
        rho = (f_x - f_new)/(0.5*d'*(lambda*d - g_x));

        % Check if step is acceptable
        if rho > 0
          % Store current iterate values
          obj.x = x_new;
          r_x = r_new;
          J_x = J_new;
          f_x = f_new;

          % Compute gradient
          g_x = J_x'*r_x;

          % Check if norm of gradient is below tolerance
          if norm(g_x, Inf) <= obj.gTolerance
            % Set status and exit loop
            obj.status = OptimizationStatus.G_TOLERANCE; break;
          end % if
          
          % Check if norm of gradient is below tolerance
          if isnan(norm(g_x, Inf))
            % Set status and exit function
            obj.status = OptimizationStatus.SINGULAR; break;
          end % if

          % Compute Gauss-Newton approximation of Hessian
          H_x = J_x'*J_x;

          % Adjust Levenberg-Marquardt parameters (decrease damping factor)
          lambda = lambda*max(1/3, 1 - (2*rho - 1)^3);
          nu = 2;
        else
          % Adjust Levenberg-Marquardt parameters (increase damping factor)
          lambda = lambda*nu;
          nu = 2*nu;
        end % if
      end % while
    end % solve
  end % methods

  % GET/SET METHODS =======================================================
  methods
    % Status --------------------------------------------------------------
    function status = getStatus(obj)
    %GETSTATUS Get solver status.
      status = obj.status;
    end % getStatus
    
    % Initial guess and solution ------------------------------------------
    function x = getSolution(obj)
    %GETSOLUTION Get solution vector.
      x = obj.x;
    end % getSolution
    
    function setInitialGuess(obj, x)
    %SETINITIALGUESS Set initial guess vector.
     
      % Validate initial guess
      validateattributes(x, ...
        {'double'}, ...
        {'real', 'size', [obj.problem.getNumberOfVariables, 1]}, ...
        'NonlinearLeastSquaresSolver.setInitialGuess', ...
        'x');
      
      % Set property
      obj.x = x;
    end % setInitialGuess
    
    % TODO - getActiveBounds
%       % Find active bounds
%       tolerance_b = 1e-4;
%       activeBounds = (obj.x >= ub - tolerance_b) - (obj.x <= lb + tolerance_b);

    function setZeroInitialGuess(obj)
    %SETZEROINITIALGUESS Set zero initial guess vector.
      obj.x = zeros([obj.problem.getNumberOfVariables, 1]);
    end % setZeroInitialGuess
    
    % Levenberg-Marquardt parameter ---------------------------------------
    function setLevenbergMarquardtParameter(obj, levenbergMarquardtParameter)
    %SETLEVENBERGMARQUARDTPARAMETER Set Levenberg-Marquardt parameter.
     
      % Validate Levenberg-Marquardt parameter
      validateattributes(levenbergMarquardtParameter, ...
        {'double'}, ...
        {'real', 'scalar', '>=', 0, '<=', 1}, ...
        'NonlinearLeastSquaresSolver.setLevenbergMarquardtParameter', ...
        'tau');
      
      % Set property
      obj.levenbergMarquardtParameter = levenbergMarquardtParameter;
    end % setLevenbergMarquardtParameter
    
    % Iteration limit -----------------------------------------------------
    function setIterationLimit(obj, iterationLimit)
    %SETITERATIONLIMIT Set solver iteration limit.
     
      % Validate property attributes
      validateattributes(iterationLimit, ...
        {'double'}, ...
        {'real', 'scalar', 'positive', 'integer'}, ...
        'NonlinearLeastSquaresSolver.setIterationLimit', ...
        'iterationLimit');
      
      % Set property
      obj.iterationLimit = iterationLimit;
    end % setIterationLimit
    
    % Termination tolerances ----------------------------------------------
    function setXTolerance(obj, xTolerance)
    %SETXTOLERANCE Set solver change in iterate tolerance.
     
      % Validate change in iterate tolerance
      validateattributes(xTolerance, ...
        {'double'}, ...
        {'real', 'scalar', 'positive'}, ...
        'NonlinearLeastSquaresSolver.setXTolerance', ...
        'xTolerance');
      
      % Set property
      obj.xTolerance = xTolerance;
    end % setXTolerance
    
    function setGTolerance(obj, gTolerance)
    %SETGTOLERANCE Set solver gradient norm tolerance.
     
      % Validate gradient tolerance
      validateattributes(gTolerance, ...
        {'double'}, ...
        {'real', 'scalar', 'positive'}, ...
        'NonlinearLeastSquaresSolver.setGTolerance', ...
        'gTolerance');
      
      % Set property
      obj.gTolerance = gTolerance;
    end % setGTolerance
  end % methods
  
  % PRIVATE STATIC METHODS ================================================
  methods (Access = private, Static = true)
    function [x, r_x, J_x] = projectOntoBounds(r, x, lb, ub)
    %PROJECTONTOBOUNDS Project residual and Jacobian onto bounds.
      
      % Loop through state vector
      for i = 1:numel(x)
        % Check if lower bound is violated
        if x(i) < lb(i)
          % Project state onto lower bound
          x(i) = lb(i);
        % Check if upper bound is violated
        elseif x(i) > ub(i)
          % Project state onto upper bound
          x(i) = ub(i);
        end % if
      end % for

      % Compute residual and Jacobian at projected iterate
      [r_x, J_x] = r.evaluate(x);
      
      % Loop through Jacobian columns
      for c = 1:size(J_x, 2)
        % Check if lower bound was violated
        if x(c) == lb(c)
          % Loop through Jacobian rows
          for r = 1:size(J_x, 1)
            % Check if Jacobian element needs projection
            if J_x(r,c) < 0
              % Project Jacobian onto lower bound
              J_x(r,c) = 0;
            end % if
          end % for
        % Check if upper bound was violated
        elseif x(c) == ub(c)
          % Loop through Jacobian rows
          for r = 1:size(J_x, 1)
            % Check if Jacobian element needs projection
            if J_x(r,c) > 0
              % Project Jacobian onto upper bound
              J_x(r,c) = 0;
            end % if
          end % for
        end % if
      end % for
    end % projectOntoBounds
  end % methods
end % classdef
