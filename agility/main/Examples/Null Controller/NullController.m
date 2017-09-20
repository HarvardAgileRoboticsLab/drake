%NULLCONTROLLER Null controller.
%
% Description:
%   This controller does nothing but can be used as a template to start
%   construction of new controllers.

% Copyright 2017 Agility Robotics
% Author: Mikhail Jones

classdef NullController < matlab.System & matlab.system.mixin.Propagates %#codegen

  % PRIVATE PROPERTIES ====================================================
  properties (Access = private)
    % Cassie inputs data structure
    inputs
    % Cassie outputs data structure
    outputs
  end % properties
  
  % PROTECTED METHODS =====================================================
  methods (Access = protected)
    % SYSTEM CLASS METHODS ================================================
    function setupImpl(obj)
    %SETUPIMPL Initialize System object.

      % Initialize the Cassie input/output data structures
      obj.inputs = CassieInputs;
      obj.outputs = CassieOutputs;
    end % setupImpl

    function inputs = stepImpl(obj, outputs)
    %STEPIMPL System output and state update equations.

      % Initialize --------------------------------------------------------
      
      % Reset the desired motor torques to zero in case they aren't defined
      obj.inputs.setZeroMotorTorques;
      
      % Update the Cassie outputs data structure
      obj.outputs.setStructure(outputs);
      
      
      % Return ------------------------------------------------------------
      
      % Return the updated Cassie inputs data structure
      inputs = obj.inputs.getStructure;
    end % stepImpl

    function resetImpl(~)
    %RESETIMPL Reset System object states.
    end % resetImpl
    
    % PROPAGATES CLASS METHODS ============================================
    function sz_1 = getOutputSizeImpl(~)
    %GETOUTPUTSIZEIMPL Get sizes of output ports.
      sz_1 = [1, 1];
    end % getOutputSizeImpl
    
    function dt_1 = getOutputDataTypeImpl(~)
    %GETOUTPUTDATATYPEIMPL Get data types of output ports.
      dt_1 = 'cassieInputsBus';
    end % getOutputDataTypeImpl
    
    function cp_1 = isOutputComplexImpl(~)
    %ISOUTPUTCOMPLEXIMPL Complexity of output ports.
      cp_1 = false;
    end % isOutputComplexImpl
    
    function flag_1 = isOutputFixedSizeImpl(~)
    %ISOUTPUTFIXEDSIZEIMPL Fixed-size or variable-size output ports.
      flag_1 = true;
    end % isOutputFixedSizeImpl
  end % methods
end % classdef