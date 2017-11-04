classdef PZTBender < DrakeSystem
    
    
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
        bact = 0.02;
        
    end
    
    methods
        
        function obj = PZTBender(name, dp, orien, ap)
            
            input_frame = MultiCoordinateFrame({CoordinateFrame('DriveVoltage', 1, {}, {'Vd'}), ...
                CoordinateFrame('ActuatorDeflections', 1, {}, {'qact'})}, [1;2]);
            output_frame = CoordinateFrame('ActuatorForce', 1, {}, {'Fact'});
            
            obj = obj@DrakeSystem(...
                0, ... % number of continuous states
                0, ... % number of discrete states
                2, ... % number of inputs: actuator voltage + deflection + def rate
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
                obj.dp.Vb = 225;
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
                obj.ap.tact = 335e-3;          %m (actuator thickness)
                obj.ap.tpzt = 135e-3;          %m (pzt thickness)
                obj.ap.wn = 3.5;            %m (nominal width)
                obj.ap.wr = 1.5;               %m (width ratio)
                obj.ap.lext = 1;            %m (extension length)
                obj.ap.lact = 9.2;          %m (actuator length)
                
                % PZT internal Stress as a function of field (see N.T. Jafferis et al. 2015)
                obj.ap.f31_eff_b = 24e-3;         %Pa*m/V (in blocked case)
                obj.ap.f31_eff_f = 32e-3;         %Pa*m/V (in free case)
                
                % PZT Young's modulus as a function of strain (see N.T. Jafferis et al. 2015)
                %                 obj.ap.Ea
                obj.ap.Emin = 38.5e3;          %sigmoid min (Pa)
                obj.ap.Emax = 81e3;            %sigmoid max (Pa)
                obj.ap.Eave = 40e3;            % average modulus (Pa)
                obj.ap.A0 = -0.3496;           %constant on cosine fit for Eave (1/m)
                obj.ap.A1 = 0.3291;            %multiplier on cosine fit for Eave (1/m)
                obj.ap.omg = 2.047;             %freq on cosine fit for Eave (unclear units)
                
                % Material properties (see N.T. Jafferis et al. 2015)
                obj.ap.Ecf = 340e3;             % Pa (CF young's modulus)
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
            df = 0.25*(1 + 2*obj.lr)*(df_num/df_den);
            obj.df = df;
            
            obj.kact = 3*(obj.ap.f31_eff_b/obj.ap.f31_eff_f)*(obj.ap.wn/obj.ap.lact^3)*obj.GF*...
                (df_den/(1+2*obj.ap.lext/obj.ap.lact));
            
        end
        
        function obj = setCFThickness(obj, tcf)
            obj.ap.tact = (obj.ap.tact - obj.tcf) + tcf;
            obj.tcf = tcf;
            
            df_den = (1/3)*obj.ap.Eave*obj.ap.tpzt*(1.5*obj.tcf^2 + 3*obj.tcf*obj.ap.tpzt + ...
                2*obj.ap.tpzt^2) + obj.ap.Ecf*obj.tcf^3/12;
            df_num = obj.ap.f31_eff_f*obj.ap.tpzt*obj.ap.lact^2*obj.eps*(obj.ap.tpzt+obj.tcf);
            obj.df = 0.25*(1 + 2*obj.lr)*(df_num/df_den);
            
            obj.kact = 3*(obj.ap.f31_eff_b/obj.ap.f31_eff_f)*(obj.ap.wn/obj.ap.lact^3)*obj.GF*...
                (df_den/(1+2*obj.ap.lext/obj.ap.lact));
            
        end
        
        
        function u0 = getDefaultInput(obj)
            u0 = [0.5*(obj.Vb - obj.Vg); 0];
        end
        
        function [u, du] = voltageToForce(obj, q, V)
            Vt = obj.dp.Vb - V;        % voltage on top plate
            Vb = V - obj.dp.Vg;        % voltage on bottom plate
            
            Ct = (0.75/obj.ap.lact)*(obj.ap.f31_eff_b + (-obj.orien*q/obj.df)*...
                (obj.ap.f31_eff_f - obj.ap.f31_eff_b))*obj.ap.wn*(obj.ap.tpzt+obj.tcf)*obj.GF;
            Cb = (0.75/obj.ap.lact)*(obj.ap.f31_eff_b + (obj.orien*q/obj.df)*...
                (obj.ap.f31_eff_f - obj.ap.f31_eff_b))*obj.ap.wn*(obj.ap.tpzt+obj.tcf)*obj.GF;
            
            dCt_dq = (0.75/obj.ap.lact)*((-obj.orien/obj.df)*(obj.ap.f31_eff_f - obj.ap.f31_eff_b))*obj.ap.wn*(obj.ap.tpzt+obj.tcf)*obj.GF;
            dCb_dq = (0.75/obj.ap.lact)*((obj.orien/obj.df)*(obj.ap.f31_eff_f - obj.ap.f31_eff_b))*obj.ap.wn*(obj.ap.tpzt+obj.tcf)*obj.GF;
            
            u = obj.orien*(Ct*Vt - Cb*Vb);            
            du = [obj.orien*(dCt_dq*Vt - dCb_dq*Vb), -obj.orien*(Ct + Cb)];
            
            
        end
        
        function [y, dy] = output(obj, t, x, u)            
            
            V = u(1);
            q = u(2);
            
            [u, du] = obj.voltageToForce(q, V);
            y = u;
            dy = [zeros(1, numel(t)), zeros(1, numel(x)), du];
        end
        
    end
    
end
