classdef ControlledPendulumFoamyPlant < FoamyPendulumPlant
    properties
        controller
    end
    
    methods
        function obj = ControlledPendulumFoamyPlant(controller, m_p)
            obj = obj@FoamyPendulumPlant(m_p);
            obj.controller = controller;
        end
        
        function [xdot, dxdot] = dynamics(obj, t, x, u)
            u = obj.controller.output(t, u, x(1:13));
            
            [xdot, dxdot] = dynamics@FoamyPlant(obj, t, x, u);
        end
        
        function setStartTime(obj, t)
            obj.controller = obj.controller.setStartTime(t);
            fprintf('Start time reset to %.5f.\n', t);
        end
    end
end