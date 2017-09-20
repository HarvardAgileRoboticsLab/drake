%NONLINEARLEASTSQUARESPROBLEM Nonlinear least squares problem object.
%
% Syntax:
%   problem = NonlinearLeastSquaresProblem(r, lb, ub);
%
% Description:
%   Least squares is the problem of finding a vector x that minimizes the
%   sum of the squares of the error. Typically this is used for regression
%   analysis or for overdetermined systems. Given a vector function r(x),
%   we want to minimize ||r(x)||, or equivalently to find
%
%     x* = argmin_x{f(x)}, where f(x) = 1/2*r(x)'*r(x)
%
%   subject to:
%
%     lb <= x <= ub
%
% TODO:
%   - Need to add upper and lower bound checks when using set methods so
%   lower and upper bounds stay consistent (lb <= ub)
%   - eorganize constructor to use set methods

% Copyright 2016-2017 Mikhail S. Jones

classdef NonlinearLeastSquaresProblem < handle %#codegen

  % PROTECTED PROPERTIES ==================================================
  properties %(SetAccess = protected)
    % Nonlinear residual function
    residualFunction
    % Decision variable lower bounds
    lowerBound
    % Decision variable upper bounds
    upperBound
  end % properties

  % PUBLIC METHODS ========================================================
  methods (Access = public)
    function obj = NonlinearLeastSquaresProblem(r, lb, ub)
    %NONLINEARLEASTSQUARESPROBLEM Construct nonlinear least squares problem.

      % Validate residual function
      validateattributes(r, ...
        {'Function'}, ...
        {'scalar'}, ...
        'NonlinearLeastSquaresProblem', ...
        'residualFunction');
      
      % Check if any lower bound constraints exist
      if isempty(lb)
        % Construct unbounded lower bound
        obj.lowerBound = -Inf(r.getNumberOfInputs, 1);
      else
        % Validate lower bound vector
        validateattributes(lb, ...
          {'double'}, ...
          {'real', 'size', [r.getNumberOfInputs, 1], '<', Inf}, ...
          'NonlinearLeastSquaresProblem', ...
          'lowerBound');
        
        % Set lower bound vector
        obj.lowerBound = lb;
      end % if
      
      % Check if any upper bound constraints exist
      if isempty(ub)
        % Construct unbounded upper bound
        obj.upperBound = Inf(r.getNumberOfInputs, 1);
      else
        % Validate upper bound vector
        validateattributes(ub, ...
          {'double'}, ...
          {'real', 'size', [r.getNumberOfInputs, 1], '>', -Inf}, ...
          'NonlinearLeastSquaresProblem', ...
          'upperBound');
        
        % Set upper bound vector
        obj.upperBound = ub;
      end % if

      % Validate lower bounds are consistent with upper bounds
      assert(...
        all(obj.lowerBound <= obj.upperBound), ...
        'Expected lowerBound to be a vector with all of the values <= upperBound.')
      
      % Set residual function
      obj.residualFunction = r;
    end % NonlinearLeastSquaresProblem

    function [r, lb, ub] = unpack(obj)
    %UNPACK Return broken apart problem structure for use in solver.

      % Unpack residual function
      r = obj.residualFunction;

      % Unpack bound constraints
      lb = obj.lowerBound;
      ub = obj.upperBound;
    end % unpack
  end % methods
  
  % GET/SET METHODS =======================================================
  methods
    function n = getNumberOfVariables(obj)
    %GETNUMBEROFVARIABLES Get the number of variables.
      n = obj.residualFunction.getNumberOfInputs;
    end % getNumberOfVariables
    
    function setLowerBound(obj, lowerBound)
    %SETLOWERBOUND Set problem lower bounds.
    
      % Validate input argument attributes
      validateattributes(lowerBound, ...
        {'double'}, ...
        {'real', 'size', [obj.residualFunction.getNumberOfInputs, 1], '<', Inf}, ...
        'NonlinearLeastSquaresProblem.setLowerBound', ...
        'lowerBound');

      % Set object properties
      obj.lowerBound = lowerBound;
    end % setLowerBound
    
    function setUpperBound(obj, upperBound)
    %SETUPPERBOUND Set problem upper bounds.
    
      % Validate input argument attributes
      validateattributes(upperBound, ...
        {'double'}, ...
        {'real', 'size', [obj.residualFunction.getNumberOfInputs, 1], '>', -Inf}, ...
        'NonlinearLeastSquaresProblem.setUpperBound', ...
        'upperBound');

      % Set object properties
      obj.upperBound = upperBound;
    end % setUpperBound
  end % methods
end % classdef
