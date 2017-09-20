%FILTER Digital filter object.
%
% Description:
%   This filter implementation uses the "Direct Form II Transposed"
%   difference equation:
%
%     a(1)*y(n) = b(1)*x(n) +  ...  + b(nb+1)*x(n-nb) ...
%               - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
%
% Notes:
%   - Coefficients are automatically normalized so a(1) is equal to one.

% Copyright 2017 Mikhail S. Jones

classdef Filter < handle %#codegen
  
  % PROTECTED PROPERTIES ==================================================
  properties (SetAccess = protected)
    % Coefficients for filtered signal
    a
    % Coefficients for unfiltered signal
    b
    % Unfiltered signal
    x
    % Filtered signal
    y
  end % properties
  
  % PUBLIC METHODS ========================================================
  methods (Access = public)
    function obj = Filter(b, a, numberOfSignals)
    %FILTER Construct digital filter object.
    
      % Set code generation settings
      coder.allowpcode('plain');
      
      % Validate input arguments
      validateattributes(b, ...
        {'numeric'}, ...
        {'vector', 'finite'}, ...
        'Filter', ...
        'b');
      validateattributes(a, ...
        {'numeric'}, ...
        {'vector', 'finite'}, ...
        'Filter', ...
        'a');
      
      % Check if number of signals was specified
      if nargin == 2
        % Set default number of signals to one
        numberOfSignals = 1;
      end % if
      
      % Set filter coefficients
      obj.b = b/a(1);
      obj.a = a/a(1);
      
      % Initialize signals
      obj.x = zeros(numberOfSignals, length(obj.b));
      obj.y = zeros(numberOfSignals, length(obj.a));
    end % Filter
    
    function updateValue(obj, value)
    %UPDATEVALUE Update signal value with digital filter.
      
      % Validate input arguments
      validateattributes(value, ...
        {'numeric'}, ...
        {'vector', 'finite', 'numel', obj.getNumberOfSignals}, ...
        'Filter.updateValue', ...
        'value');
      
      % Shift unfiltered signal array
      obj.x = circshift(obj.x, 1, 2);
      
      % Shift filtered signal array
      obj.y = circshift(obj.y, 1, 2);
      
      % Update value of unfiltered signal at current point
      obj.x(:,1) = value;
      
      % Reshape coefficients
      A = repmat(obj.a, obj.getNumberOfSignals, 1);
      B = repmat(obj.b, obj.getNumberOfSignals, 1);
      
      % Update value of unfiltered signal at current point
      if length(obj.a) == 1
        obj.y(:,1) = sum(B.*obj.x, 2);
      else
        obj.y(:,1) = sum(B.*obj.x, 2) - sum(A(:,2:end).*obj.y(:,2:end), 2);
      end % if
    end % updateValue
  end % methods
  
  % GET/SET METHODS =======================================================
  methods
    % Number of signals ---------------------------------------------------
    function numberOfSignals = getNumberOfSignals(obj)
    %GETNUMBEROFSIGNALS Get the number of signals.
      numberOfSignals = size(obj.x, 1);
    end % getNumberOfSignals
    
    % Signal value --------------------------------------------------------
    function value = getValue(obj)
    %GETVALUE Get signal value.
      value = obj.y(:,1);
    end % getValue
    
    function setUnfilteredValue(obj, value)
    %SETUNFILTEREDVALUE Set unfiltered signal value.
    
      % Validate input arguments
      validateattributes(value, ...
        {'numeric'}, ...
        {'vector', 'finite', 'numel', obj.getNumberOfSignals}, ...
        'Filter.setUnfilteredValue', ...
        'value');
      
      % Set unfiltered signal history
      obj.x = repmat(value, 1, length(obj.b));
    end % setUnfilteredValue
    
    function setFilteredValue(obj, value)
    %SETFILTEREDVALUE Set filtered signal value.
    
      % Validate input arguments
      validateattributes(value, ...
        {'numeric'}, ...
        {'vector', 'finite', 'numel', obj.getNumberOfSignals}, ...
        'Filter.setFilteredValue', ...
        'value');
      
      % Set filtered signal history
      obj.y = repmat(value, 1, length(obj.a));
    end % setFilteredValue
  end % methods
  
  % STATIC METHODS ========================================================
  methods (Static = true)
    function C = cascade(A, B)
    %CASCADE Cascade filter structure.
    
      % Validate input arguments
      validateattributes(A, ...
        {'Filter'}, ...
        {'scalar'}, ...
        'Filter.cascade', ...
        'A');
      validateattributes(B, ...
        {'Filter'}, ...
        {'scalar'}, ...
        'Filter.cascade', ...
        'B');
      
      % Check that filter have same number of signals
      assert(A.getNumberOfSignals == B.getNumberOfSignals, ...
        'Filters must have the same number of signals.');
      
      % Compute cascaded coefficients
      a = conv(A.a, B.a);
      b = conv(A.b, B.b);
      
      % Construct cascaded filter
      C = Filter(b, a, A.getNumberOfSignals);
    end % cascade
    
    function C = parallel(A, B)
    %PARALLEL Parallel filter structure.
      
      % Validate input arguments
      validateattributes(A, ...
        {'Filter'}, ...
        {'scalar'}, ...
        'Filter.parallel', ...
        'A');
      validateattributes(B, ...
        {'Filter'}, ...
        {'scalar'}, ...
        'Filter.parallel', ...
        'B');
      
      % Check that filter have same number of signals
      assert(A.getNumberOfSignals == B.getNumberOfSignals, ...
        'Filters must have the same number of signals.');
      
      % Compute maximum size of each coefficient vector
      na = max(numel(A.a), numel(B.a));
      nb = max(numel(A.b), numel(B.b));
      
      % Initialize new coefficient vectors
      a = zeros(1,na);
      b = zeros(1,nb);
      
      % Compute parallel coefficients
      a(1:numel(A.a)) = a(1:numel(A.a)) + A.a;
      a(1:numel(B.a)) = a(1:numel(B.a)) + B.a;
      b(1:numel(A.b)) = b(1:numel(A.b)) + A.b;
      b(1:numel(B.b)) = b(1:numel(B.b)) + B.b;
      
      % Construct parallel filter
      C = Filter(b, a, A.getNumberOfSignals);
    end % parallel
    
    function B = copy(A)
    %COPY Copy filter structure to new object.
    
      % Validate input arguments
      validateattributes(A, ...
        {'Filter'}, ...
        {'scalar'}, ...
        'Filter.copy', ...
        'A');
      
      % Construct copied filter
      B = Filter(A.b, A.a, A.getNumberOfSignals);
    end % copy
  end % methods
end % classdef