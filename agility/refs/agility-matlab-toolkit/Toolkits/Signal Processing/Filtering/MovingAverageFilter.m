%MOVINGAVERAGEFILTER Moving average filter object.

% Copyright 2017 Mikhail S. Jones

classdef MovingAverageFilter < Filter %#codegen
  
  % PUBLIC METHODS ========================================================
  methods (Access = public)
    function obj = MovingAverageFilter(windowSize, varargin)
    %MOVINGAVERAGEFILTER Construct moving average filter object.
    
      % Set code generation settings
      coder.allowpcode('plain');
      
      % Compute filter coefficients
      b = repmat(1/windowSize, 1, windowSize);
      a = 1;
      
      % Call superclass constructor
      obj = obj@Filter(b, a, varargin{:});
    end % MovingAverageFilter
  end % methods
end % classdef