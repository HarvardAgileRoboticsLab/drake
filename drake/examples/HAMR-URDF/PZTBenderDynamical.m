classdef PZTBenderDynamical < DrakeSystem
    
    
    properties
        
        name = '';            %name
        
        ap = struct();        % actuator parameters
        dp = struct();        % drive signal
        
        % orientation (1 = pointed forward, -1 = pointed backwards)
        orien
        
        % Derived quantities
        tcf;        % carbon fiber thickness
        lr;         % length ratio
        eps;        % total field
        GF;         % geometry factor
        df;         % free deflection
        kact;       % stiffness
        bact = 0.01; %0.00263; % damping (TODO: GET GOOD VALUE)
        
        % Hamr deflection joints
        %         obj.
        
    end
    
    methods
        
        function obj = PZTBender(name, dp, orien, ap)
            
            input_frame = MultiCoordinateFrame({CoordinateFrame('DriveVoltage', 1, {}, {'Vd'}), ...
                CoordinateFrame('ActuatorDeflectionandRate', 2, {})}, [1;2;2]);
            %             input_frame = CoordinateFrame('DriveVoltage', 1, {}, {'Vd'});
            output_frame = CoordinateFrame('ActuatorForce', 1, {}, {'Fact'});
            
            obj = obj@DrakeSystem(...
                0, ... % number of continuous states
                0, ... % number of discrete states
                3, ... % number of inputs: actuator voltage + deflection + def rate
                1, ... % number of outputs: Actuator force/stiffness
                true, ... % direct feedthrough
                true); % time invariant
            
            obj = obj.setInputFrame(input_frame);
            obj = obj.setOutputFrame(output_frame);
            
            if nargin < 1
                obj.name = 'PZTBender';
            else
                obj.name = name;
            end
            
            if nargin < 2
                obj.dp.Vb = 200;
                obj.dp.Vg = 0;
            else
                obj.dp = dp;
            end
            
            if nargin < 3
                obj.orien = 1;      % default is forward
            else
                obj.orien = orien;
            end
            
            % Set limits
            %obj = obj.setInputLimits([obj.dp.Vg; -Inf], [obj.dp.Vb; Inf]);
            
            if nargin < 4 % Default geometry for HAMR
                obj.ap.tact = 325e-6;          %m (actuator thickness)
                obj.ap.tpzt = 135e-6;          %m (pzt thickness)
                obj.ap.wn = 3.5e-3;            %m (nominal width)
                obj.ap.wr = 1.5;               %m (width ratio)
                obj.ap.lext = 1e-3;            %m (extension length)
                obj.ap.lact = 9.2e-3;          %m (actuator length)
                
                % PZT internal Stress as a function of field (see N.T. Jafferis et al. 2015)
                obj.ap.f31_eff_f = 32;                     % Pa*m/Volts (sigmoid max)
                obj.ap.f31_eff_b = 24;                     % Pa*m/Volts (sigmoid min)
                
                % PZT Young's modulus as a function of strain (see N.T. Jafferis et al. 2015)
                obj.ap.Eave = 40e9;                % average modulus
                
                % Material properties (see N.T. Jafferis et al. 2015)
                obj.ap.Ecf = 340e9;             % Pa (CF young's modulus)
            else
                obj.ap = ap;
            end
            
            %  ----------- fillin in derived actuator properties -------------
            obj.tcf = obj.ap.tact - 2*obj.ap.tpzt;            % cf thickness
            obj.lr = obj.ap.lext/obj.ap.lact;                 % length ratio
            obj.eps = (obj.dp.Vb -obj.dp.Vg)/obj.ap.tpzt;        % total field
            
            % geometric factor
            GF_num = 8*(1-obj.ap.wr)^3*(1 + 2*obj.lr);
            GF_den = -6*(obj.ap.wr - 1)*(-3 + 4*obj.lr*(obj.ap.wr - 1)+2*obj.ap.wr) ...
                + 3*(-2 + 2*obj.lr*(obj.ap.wr-1)+obj.ap.wr)^2*log((2-obj.ap.wr)/obj.ap.wr);
            obj.GF = GF_num/GF_den;
            
            % free deflection
            df_den = (1/3)*obj.ap.Eave*obj.ap.tpzt*(1.5*obj.tcf^2 + 3*obj.tcf*obj.ap.tpzt + ...
                2*obj.ap.tpzt^2) + obj.ap.Ecf*obj.tcf^3/12;
            df_num = obj.ap.f31_eff_f*obj.ap.tpzt*obj.ap.lact^2*obj.eps*(obj.ap.tpzt+obj.tcf);
            obj.df = 0.25*(1 + 2*obj.lr)*(df_num/df_den);
            
            % actuator stiffnes
            obj.kact = 3.0*obj.ap.wn*df_den*obj.GF/(obj.ap.lact^3*(1+2*obj.lr));
                        %             obj.bact = 3*obj.kact;
            
        end
        
        
        function [xdot, dxdot] = dynamics(obj,t,x,u)
            
            Vt = obj.dp.Vb - u(1);        % voltage on top plate
            Vb = u(1) - obj.dp.Vg;        % voltage on bottom plate
            q = x(1);
            qd = x(2); 

            Ft = (0.75*Vt/obj.ap.lact)*obj.ap.f31_eff_b*obj.ap.wn*(obj.ap.tpzt+obj.tcf)*obj.GF;
            Fb = (0.75*Vb/obj.ap.lact)*obj.ap.f31_eff_b*obj.ap.wn*(obj.ap.tpzt+obj.tcf)*obj.GF;
            
            qdd = (obj.orien*(Ft - Fb) - obj.kact*q - obj.bact*qd)*(1/obj.m); 
            xdot = [qd; qdd]; 
            
        end
        
        function u0 = getDefaultInput(obj)
            u0 = [0.5*(obj.Vb - obj.Vg); 0];
        end
        
        % here's where the real stuff happens
        function [y, dy] = output(obj, t, x, u)
            y = x; 
            if nargout > 1
                dy = eye(size(x)); 
            end
        end
    end
    
end
