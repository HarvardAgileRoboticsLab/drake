classdef PublishStateLCM < matlab.System & matlab.system.mixin.Propagates %#codegen

  % PRIVATE PROPERTIES ====================================================
  properties (Access = private)
    outputs
  end % properties
  
  % PROTECTED METHODS =====================================================
  methods (Access = protected)
    % SYSTEM CLASS METHODS ================================================
    function setupImpl(obj)
    %SETUPIMPL Initialize System object.

      % Initialize the Cassie input/output data structures
      obj.outputs = CassieOutputs;
    end % setupImpl

    function outputs = stepImpl(obj, outputs)
    %STEPIMPL System output and state update equations.

      coder.extrinsic('lcmPublishHelper') 
      lcmPublishHelper(outputs);
      
      % Update the Cassie outputs data structure
      obj.outputs.setStructure(outputs);
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
      dt_1 = 'cassieOutputsBus';
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