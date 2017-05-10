classdef OLControl < DrakeSystem
    
    
    properties
        
    end
    
    methods
        
        function obj = OLControl()
            
            %input_frame = MultiCoordinateFrame({CoordinateFrame('DriveVoltage', 1, {}, {'Vd'}), ...
            %    CoordinateFrame('ActuatorDeflection', 1, {})}, [1;2]);
%             output_frame = CoordinateFrame('DriveVoltage', 1, {}, {'Vd'});
%             output_frame = CoordinateFrame('ActuatorForce', 1, {}, {'Fact'});
            
            obj = obj@DrakeSystem(...
                1, ... % number of continuous states: dummy state
                0, ... % number of discrete states
                1, ... % number of inputs: 1
                1, ... % number of outputs: drive voltage
                true, ... % direct feedthrough
                false); % time invariant
        end
           
           
        function xdot = dynamics(obj,t,x,u)
            
            xdot = 75*2*pi*cos(2*pi*t); 
                        %dxdot = [1, zeros(1, 1 + numel(x) + numel(u))];
            
        end
        
        % here's where the real stuff happens
        function y = output(obj, t, x, u)            
            tramp = 3;
            
            if t <= tramp
                ramp = t/tramp;
            else 
                ramp = 1;
            end
            
            y = ramp*x + 75;
            fprintf('Vact should be: %f\n', y);  
            fprintf('And it is %f\n', u); 
            
        end

    end
    
end
