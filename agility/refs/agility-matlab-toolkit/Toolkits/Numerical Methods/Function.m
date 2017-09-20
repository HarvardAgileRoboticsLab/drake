%FUNCTION Abstract function object.
%
% Description:
%   This implements an abstract class for building functionality on top
%   of functions.

% Copyright 2016 Mikhail S. Jones

classdef (Abstract = true) Function < handle %#codegen

  % PROTECTED PROPERTIES ==================================================
  properties (Access = protected)
    % Function name
    name
    % Number of input variables
    numberOfInputs
    % Number of output variables
    numberOfOutputs
  end % properties

  % ABSTRACT METHODS ======================================================
  methods (Abstract = true)
    % Evaluate function
    [f, J] = evaluate(obj, x)
  end % methods
  
  % PUBLIC METHODS ========================================================
  methods (Access = public)
    function obj = Function(name, numberOfInputs, numberOfOutputs)
    %FUNCTION Function object constructor.
    
      % Set code generation settings
      coder.allowpcode('plain');
      
      % Validate function name
      validateattributes(name, ...
        {'char'}, ...
        {'vector'}, ...
        'Function', ...
        'name');
      
      % Set function name
      obj.name = name;
 
      % Validate number of inputs
      validateattributes(numberOfInputs, ...
        {'double'}, ...
        {'real', 'scalar', 'integer', 'positive'}, ...
        'Function', ...
        'numberOfInputs');
      
      % Set number of inputs
      obj.numberOfInputs = zeros(numberOfInputs, 1);
      
      % Validate number of outputs
      validateattributes(numberOfOutputs, ...
        {'double'}, ...
        {'real', 'scalar', 'integer', 'positive'}, ...
        'Function', ...
        'numberOfOutputs');
      
      % Set number of outputs
      obj.numberOfOutputs = zeros(numberOfOutputs, 1);
    end % Function
  end % methods
  
  % GET/SET METHODS =======================================================
  methods
    % Input/Output Size ---------------------------------------------------
    function n = getNumberOfInputs(obj)
    %GETNUMBEROFINPUTS Get the number of inputs. 
      n = numel(obj.numberOfInputs);
    end % getNumberOfInputs
    
    function n = getNumberOfOutputs(obj)
    %GETNUMBEROFOUTPUTS Get the number of outputs. 
      n = numel(obj.numberOfOutputs);
    end % getNumberOfOutputs
    
    % Name ----------------------------------------------------------------
    function name = getName(obj)
    %GETNAME Get name of this function.
      name = obj.name;
    end % getName
    
    function setName(obj, name)
    %SETNAME Set name of this function.
    
      % Validate name attributes
      validateattributes(name, ...
        {'char'}, ...
        {'vector'}, ...
        'Function.setName', ...
        'name');
      
      % Set name property
      obj.name = name;
    end % setName
  end % methods
end % classdef
