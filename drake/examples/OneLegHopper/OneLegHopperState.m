classdef OneLegHopperState < SingletonCoordinateFrame
  methods
    function obj=OneLegHopperState(r)
      obj = obj@SingletonCoordinateFrame('OneLegHopperState',r.getNumStates,'x');
    end
  end
end
