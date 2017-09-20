%BACKWARDDIFFERENCINGFILTER Backward differencing filter object.

% Copyright 2017 Mikhail S. Jones

classdef BackwardDifferencingFilter < Filter %#codegen
  
  % PUBLIC METHODS ========================================================
  methods (Access = public)
    function obj = BackwardDifferencingFilter(numberOfPoints, sampleTime, varargin)
    %BACKWARDDIFFERENCINGFILTER Construct backward differencing filter object.
    
      % Set code generation settings
      coder.allowpcode('plain');
      
      % Validate input arguments
      validateattributes(numberOfPoints, ...
        {'numeric'}, ...
        {'scalar', 'integer', '>=', 2}, ...
        'BackwardDifferencingFilter', ...
        'numberOfPoints');
      validateattributes(sampleTime, ...
        {'numeric'}, ...
        {'scalar', 'positive'}, ...
        'BackwardDifferencingFilter', ...
        'sampleTime');
      
      % Initialize coefficient array
      % b = zeros(1, numberOfPoints);
      assert(numberOfPoints == 2, ...
        'More than two points not supported yet.');
            
      % Call superclass constructor
      obj = obj@Filter([1 -1], sampleTime, varargin{:});
    end % BackwardDifferencingFilter
  end % methods
end % classdef