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
            rmax0 = 0.67*(rp-rh);
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
            
            %Smooth fit to data in Khan and Nahon (2015)
            function vs = vfit(pos)
                x = (-pos(1)-x0)/D0;
                r = norm(pos(2:3));
                
                vs = vmax(x)*exp(-(r-rmax(x))/(a3(x)*rmax0 + b3(x)*(x*D0-R0)));
            end
            
            function vm = vmax(x)
                a = 0.3197;
                b = -0.527;
                c = -1.413;
                d =  0.9509;
                
                vm = v0*a*tanh(b*x-c)+d;
            end
            
            function rm = rmax(x)
                a = -0.5138;
                b = 0.7461;
                c = 1.956;
                d = 0.4809;
                
                rm = rmax0*a*tanh(b*x-c)+d;
            end
            
            function a = a3(x)
                a = -0.1546*x + 0.9285;
            end
            
            function b = b3(x)
                b = -0.002997*x^2 + 0.03692*x + 0.1298;
            end
            
        end
        
        
        
    end
    
end

