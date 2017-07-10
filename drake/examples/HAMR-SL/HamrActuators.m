classdef HamrActuators < DrakeSystem
    
    
    properties
        names
        dummy_bender = [];
        nact;
    end
    
    methods
        
        function obj = HamrActuators(nact, names, orien, dp)     
        % TODO: add option to change actuator shape by passing gp
            
            if numel(names) ~= nact
                error('Need a name for each acutator mofo')
            end            
            
            if nargin < 3
                orien = ones(nact, 1); 
            else
                if numel(orien) ~= nact
                    disp('Orien is wrong length, setting to ones')            
                    orien = ones(nact,1);                     
                end
            end
            
            if nargin < 4
                dp.Vb = 200;
                dp.Vg = 0;                                       
            end

            input_frame = MultiCoordinateFrame({CoordinateFrame('DriveVoltage', nact, {}, names), ...
                CoordinateFrame('ActuatorDeflectionandRate', 2*nact)} ...
                ,[ones(nact,1); 2*ones(2*nact,1)]);
            
            output_frame = CoordinateFrame('ActuatorForce', nact, {}, names);
            
            obj = obj@DrakeSystem( ...
                0, ... % number of continuous states: dummy state
                0, ... % number of discrete states
                3*nact, ... % number of inputs: actuator voltage + deflection and rate
                nact, ... % number of outputs: Actuator force/stiffness
                true, ... % direct feedthrough
                true); % time invariant
            
            obj = obj.setInputFrame(input_frame);
            obj = obj.setOutputFrame(output_frame);
            
            obj.names = names;
            for i = 1:nact
                obj.dummy_bender = [obj.dummy_bender, PZTBender(names{i}, dp,orien(i))];
            end
            obj.nact = nact;            
        end
        
        % here's where the real stuff happens
        function y = output(obj, t, x, u)
            y = zeros(obj.nact, 1);
            for i = 1:obj.nact
                y(i) = obj.dummy_bender(i).output(t, x, [u(i);...
                    u(i+obj.nact); u(i+2*obj.nact)]);
            end
        end
        
        function u0 = getDefaultInput(obj)
            u0 = [0.5*(obj.Vb - obj.Vg)*ones(obj.nact, 1); 
                zeros(obj.nact,1)]; 
        end
        
    end
    
end
