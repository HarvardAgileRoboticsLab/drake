classdef McFoamy < Airplane
    %MCFOAMY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        m = 0.484;
        J = [0.003922 0.000303 0.000441;
             0.000303 0.015940 -0.000030;
             0.000441 -0.000030 0.019340];
        geom; %Holds a geometry struct that describes the airplane
    end
    
    methods
        function obj = McFoamy()
            obj = obj@Airplane();
            obj.geom = obj.geomSetup();
        end
        
        function [F,T,dF,dT] = aerodynamics(obj,t,q,v,w,thr,ail,ele,rud)
            F = [0 0 0]';
            T = [0 0 0]';
        end
        
        function g = geomSetup(obj)
            g = struct();
            
            g.prop = struct();
            g.wing = struct();
            g.hTail = struct();
            g.vTail = struct();
            g.body = struct();
            
            g.cg = [-270 0 5.89]'; %position of CG relative to propeller
            
            g.prop.d = 254.0; %overall prop diameter (blade + hub)
            g.prop.dh = 12.0; %diameter of prop hub
            
            g.wing.span = [57.26 58.72 59.29 62.58 62.58 63.25 61.90]';
            g.wing.area = [14499.66 13977.71 13198.25 12974.40 11992.21 11122.51 9913.9]';
            g.wing.chord = [253.4000  237.5800  222.6900  207.4200  191.7400  175.9700  160.2800]';
            g.wing.flap = [0  102.0400   99.2300   96.3600   93.4100   90.4400   87.4900]';
            g.wing.pos = {[-217.0700   33.3800         0]', ... 
                          [-211.3900   91.7600         0]', ...
                          [-215.4300  151.1400         0]', ...
                          [-212.1800  212.0100         0]', ...
                          [-208.8500  274.5500         0]', ...
                          [-205.5000  337.4200         0]', ...
                          [-202.2600  400.2100         0]'};
                      
            g.hTail.span = [66.5500   51.0900   62.0200]';
            g.hTail.area = [7241.76    7452.58    8892.67]';
            g.hTail.chord = [117.4100  145.6400  143.1600]';
            g.hTail.flap = [81.8700  111.4900  143.1600]';
            g.hTail.pos = {[-711.2700   39.5300         0]', ...
                           [-719.7200   93.3300         0]', ...
                           [-715.1500  151.5700         0]'};
            
            g.vTail.span = [46.3000   62.1900   69.3100   38.1000]';
            g.vTail.area = [9291.95  11997.38  11715.82   5298.95]';
            g.vTail.chord = [200.6200  192.9500  169.4900  141.1400]';
            g.vTail.flap = [120.3600  112.6900  103.4300  141.1400]';
            g.vTail.pos = {[-733.9400         0   23.5400], ...
                           [-735.8100         0  -29.9600], ...
                           [-744.1400         0  -94.6300], ...
                           [-760.7600         0 -148.5600]'};
                       
            g.body.span = [35.8700   35.8700   37.2600   37.2600]';
            g.body.area = [23018.14  23018.14  23910.11  23910.11]';
            g.body.chord = [641.7100  641.7100  641.7100  641.7100]';
            g.body.flap = [0     0     0     0]';
            g.body.pos = {[-205.0900         0   53.8100], ...
                          [-205.0900         0   17.9400], ...
                          [-205.0900         0  -18.6300], ...
                          [-205.0900         0  -55.8900]'};        
        end
        
        function vstream = slipstream(obj, vi0)
            %Implements the propeller slipstream model from Khan and Nahon (2015)
            
            rp = obj.geom.prop.d/(2e-3);
            rh = obj.geom.prop.dh/(2e-3);
            
            x0 = 1.528*rp; %eq. 11
            R0 = 0.74*rp;  %eq. 12
            D0 = 2*R0;
            Rm0 = 0.67*(rp-rh);
            v0 = 2*vi0*(1.46/1.59); %efflux velocity
            
            
            vstream.wing = zeros(length(obj.geom.wing.pos),1);
            for k = 1:length(obj.geom.wing.pos)
                vstream.wing(k) = vfit(obj.geom.wing.pos{k});
            end
            
            vstream.hTail = zeros(length(obj.geom.hTail.pos),1);
            for k = 1:length(obj.geom.hTail.pos)
                vstream.hTail(k) = vfit(obj.geom.hTail.pos{k});
            end
            
            vstream.vTail = zeros(length(obj.geom.vTail.pos),1);
            for k = 1:length(obj.geom.vTail.pos)
                vstream.vTail(k) = vfit(obj.geom.vTail.pos{k});
            end
            
            vstream.body = zeros(length(obj.geom.body.pos),1);
            for k = 1:length(obj.geom.body.pos)
                vstream.body(k) = vfit(obj.geom.body.pos{k});
            end
            
            %TODO: replace this with a smooth interpolant
            function vs = vfit(pos)
                x = -pos(1);
                if x < x0
                    vs = vi0*(1+x/rp)/sqrt(1 + (x*x/(rp*rp))); %eq. 1
                elseif x < 1.7*D0
                    r = norm(pos(2:3));
                    vmax = v0*(1.24 - 0.0765*(x-x0)/D0); %eq. 13
                    Rm = Rm0*(1 - 0.1294*(x-x0)/D0);     %eq. 14
                    vs = vmax*exp(-((r-Rm)/(0.8839*Rm0 + 0.1326*(x-x0-R0)))^2); %eq. 15
                elseif x < 4.25*D0
                    r = norm(pos(2:3));
                    vmax = v0*(1.37 - 0.1529*(x-x0)/D0);  %eq. 16
                    Rm = Rm0*(1.3 - 0.3059*(x - x0)/D0);  %eq. 17
                    vs = vmax*exp(-((r-Rm)/(0.5176*Rm0 + 0.2295*(x-x0-R0)))^2); %eq. 18
                else
                    r = norm(pos(2:3));
                    vmax = v0*(0.89 - 0.04*(x-x0)/D0);    %eq. 19
                    vs = vmax*exp(-(r/(0.2411*(x-x0)))^2);%eq. 21
                end
            end
        end
        
        
        
    end
    
end

