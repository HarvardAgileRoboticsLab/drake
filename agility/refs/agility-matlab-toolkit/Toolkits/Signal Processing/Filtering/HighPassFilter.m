%HIGHPASSFILTER High pass filter object.

% Copyright 2017 Mikhail S. Jones

classdef HighPassFilter < Filter %#codegen
  
  % PUBLIC METHODS ========================================================
  methods (Access = public)
    function obj = HighPassFilter(varargin)
    %HIGHPASSFILTER Construct high pass filter object.
    
      % Set code generation settings
      coder.allowpcode('plain');
      
      % Call superclass constructor
      obj = obj@Filter([1 -1], [1 -1], varargin{:});
    end % HighPassFilter
    
    function setSmoothingFactor(obj, alpha)
    %SETSMOOTHINGFACTOR Set smoothing factor.
      
      % Set filter coefficients
      obj.b = [alpha, -alpha];
      obj.a = [1, -alpha];
    end % setSmoothingFactor
    
    function setTimeConstant(obj, tau, dt)
    %SETTIMECONSTANT Set time constant.
      
      % Compute smoothing factor
      alpha = tau/(tau + dt);
      
      % Set filter coefficients
      obj.b = [alpha, -alpha];
      obj.a = [1, -alpha];
    end % setTimeConstant
    
    function setCutoffFrequency(obj, f_c, dt)
    %SETCUTOFFFREQUENCY Set cutoff frequency.
    
      % Compute time constant
      tau = 1/(2*pi*f_c);
      
      % Compute smoothing factor
      alpha = tau/(tau + dt);
      
      % Set filter coefficients
      obj.b = [alpha, -alpha];
      obj.a = [1, -alpha];
    end % setCutoffFrequency
  end % methods
end % classdef