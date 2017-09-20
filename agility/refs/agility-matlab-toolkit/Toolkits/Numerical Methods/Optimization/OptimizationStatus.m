%OPTIMIZATIONSTATUS Optimization solver status enumeration.
%
% Syntax:
%   status = OptimizationStatus.ENUMERATION;
%
% Description:
%   Specifies status for optimization solvers. Negative number correspond 
%   to unsuccessful results and positive numbers correspond to successful 
%   results. Use enumeration command to view available enumerations.

% Copyright 2016 Mikhail S. Jones

classdef OptimizationStatus < int8

  % ENUMERATION ===========================================================
  enumeration
    % Problem is infeasible
    INFEASIBLE (-3)
    % Singular or badly scaled matrix
    SINGULAR (-2)
    % Iteration limit was reached
    ITERATION_LIMIT (-1)
    % Iterating
    ITERATE (0)
    % Change in cost function tolerance met
    F_TOLERANCE (1)
    % Change in independent variables tolerance met
    X_TOLERANCE (2)
    % Norm of gradient tolerance met
    G_TOLERANCE (3)
    % Total relative error tolerance met
    PHI_TOLERANCE (4)
  end % enumeration

  % PUBLIC METHODS ========================================================
  methods (Access = public)
    function string = char(obj)
    %CHAR Convert enumeration into a string.
    %
    % Syntax:
    %   string = char(status);

      % Parse enumeration
      switch obj
        case OptimizationStatus.INFEASIBLE
          string = 'Problem is infeasible';
        case OptimizationStatus.SINGULAR
          string = 'Matrix is singular, results may be inaccurate.';
        case OptimizationStatus.ITERATION_LIMIT
          string = 'Iteration limit was exceeded';
        case OptimizationStatus.ITERATE
          string = 'Solver is iterating';
        case OptimizationStatus.F_TOLERANCE
          string = 'Change in cost function below specified tolerance';
        case OptimizationStatus.X_TOLERANCE
          string = 'Change in iterate below specified tolerance';
        case OptimizationStatus.G_TOLERANCE
          string = 'Norm of gradient below specified tolerance';
        case OptimizationStatus.PHI_TOLERANCE
          string = 'Total relative error below specified tolerance';
        otherwise
          error('Should never get here');
      end % switch
    end % char
  end % methods
end % classdef
