classdef HamrActuators < DrakeSystem
    
    
    properties
        names
        dummy_bender = [];
        nact;
        orien = [-1; -1; 1; 1; -1; -1; 1; 1];       % default orientation
    end
    
    methods
        
        function obj = HamrActuators(nact, names, orien, dp)
            % TODO: add option to change actuator shape by passing gp
            
            input_frame = MultiCoordinateFrame({CoordinateFrame('DriveVoltage', nact, {}, names), ...
                CoordinateFrame('ActuatorDeflection', nact)} ...
                ,[ones(nact,1); 2*ones(nact,1)]);
            
            output_frame = CoordinateFrame('ActuatorForce', nact, {}, names);
            
            obj = obj@DrakeSystem( ...
                0, ... % number of continuous states: dummy state
                0, ... % number of discrete states
                2*nact, ... % number of inputs: actuator voltage + deflection
                nact, ... % number of outputs: Actuator force/stiffness
                true, ... % direct feedthrough
                true); % time invariant
            
            obj = obj.setInputFrame(input_frame);
            obj = obj.setOutputFrame(output_frame);
            
            
            if numel(names) ~= nact
                error('Need a name for each acutator mofo')
            end
            
            if nargin < 3
                orien = ones(nact, 1);
            else
                if isempty(orien)
                    %                     disp('Orien is wrong length, setting to ones')
                    disp('Using default actuator orientation')
                    orien = obj.orien;
                end
            end
            
            if nargin < 4
                dp.Vb = 225;
                dp.Vg = 0;
            end           
               
            
            obj.names = names;
            for i = 1:nact
                obj.dummy_bender = [obj.dummy_bender, PZTBender(names{i}, dp,orien(i))];
            end
            obj.nact = nact;
            %             obj = obj.setInputLimits([(dp.Vg-1)*ones(nact, 1); -Inf(nact,1)], ...
            % [dp.Vb*ones(nact, 1); Inf(nact,1)]);
            
        end
        
        function [y, dy] = output(obj, t, x, u)            
           
            xin = [t; x; u];
            [y,dy] = output_fun(obj,xin);
%             
%             dy_fd = zeros(size(dy));
%             dxin = 1e-6*eye(length(xin));
%             for k = 1:length(xin)
%                 dy_fd(:,k) = (output_fun(obj,xin+dxin(:,k)) - ...
%                     output_fun(obj,xin-dxin(:,k)))/2e-6;
%             end
%             
%             disp('Output derivative error:');
%             disp(max(abs(dy_fd(:)-dy(:))));
        end            
        
        % here's where the real stuff happens
        function [y, dy] = output_fun(obj, xin)
            
            t = xin(1);
            x = [];
            u = xin(1+(1:2*obj.nact)); 
            
            y = zeros(obj.nact, 1);
            dy = zeros(obj.nact, 1+2*obj.nact);
            for i = 1:obj.nact
                [y(i), dyi] = obj.dummy_bender(i).output(t, x, [u(i);...
                    u(i+obj.nact)]);
                dy(i, 1+[i, i + obj.nact]) = dyi(2:end);
            end
            %             fprintf('Time: %f \n', t);
            
        end
        
        function u0 = getDefaultInput(obj)
            u0 = [0.5*(obj.Vb - obj.Vg)*ones(obj.nact, 1);
                zeros(obj.nact,1)];
        end
        
        %         function x0 = getInitialState(obj)
        %             x0 = 0;
        %         end
        
    end
    
end
