%LOWPASSFILTER Low pass filter object.

% Copyright 2017 Mikhail S. Jones

classdef LowPassFilter < Filter %#codegen
  
  % PUBLIC METHODS ========================================================
  methods (Access = public)
    function obj = LowPassFilter(varargin)
    %LOWPASSFILTER Construct low pass filter object.
    
      % Set code generation settings
      coder.allowpcode('plain');
      
      % Call superclass constructor
      obj = obj@Filter(1, [1 0], varargin{:});
    end % LowPassFilter
    
    function setSmoothingFactor(obj, alpha)
    %SETSMOOTHINGFACTOR Set smoothing factor.
      
      % Set filter coefficients
      obj.b = alpha;
      obj.a = [1, alpha - 1];
    end % setSmoothingFactor
    
    function setTimeConstant(obj, tau, dt)
    %SETTIMECONSTANT Set time constant.
      
      % Compute smoothing factor
      alpha = dt/(tau + dt);
      
      % Set filter coefficients
      obj.b = alpha;
      obj.a = [1, alpha - 1];
    end % setTimeConstant
    
    function setCutoffFrequency(obj, f_c, dt)
    %SETCUTOFFFREQUENCY Set cutoff frequency.
    
      % Compute time constant
      tau = 1/(2*pi*f_c);
      
      % Compute smoothing factor
      alpha = dt/(tau + dt);
      
      % Set filter coefficients
      obj.b = alpha;
      obj.a = [1, alpha - 1];
    end % setCutoffFrequency
  end % methods
end % classdef