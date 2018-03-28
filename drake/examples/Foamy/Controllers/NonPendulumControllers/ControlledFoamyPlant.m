classdef ControlledFoamyPlant < FoamyPlant
    properties
        controller
    end
    
    methods
        function obj = ControlledFoamyPlant(controller)
            obj = obj@FoamyPlant();
            obj.controller = controller;
        end
        
        function [xdot, dxdot] = dynamics(obj, t, x, u)
            u = obj.controller.output(t, u, x);
            
            [xdot, dxdot] = dynamics@FoamyPlant(obj, t, x, u);
        end
        
        function setStartTime(obj, t)
            obj.controller = obj.controller.setStartTime(t);
            fprintf('Start time reset to %.5f.\n', t);
        end
    end
end