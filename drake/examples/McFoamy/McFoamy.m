classdef McFoamy < Airplane
    %MCFOAMY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        m = 0.484;
        J = [0.003922 0.000303 0.000441;
             0.000303 0.015940 -0.000030;
             0.000441 -0.000030 0.019340];
        Jinv
        Cd0 = 0.02;
        Cd90 = 1.98;
        
        geom; %Holds a geometry struct that describes the airplane
    end
    
    methods
        function obj = McFoamy()
            obj = obj@Airplane();
            obj.geom = obj.geomSetup();
            obj.Jinv = inv(obj.J);
        end
        
        function [F,M,dF,dM] = aerodynamics(obj,t,q,v,w,thr,ail,ele,rud)
            
            %rotate velocity into body frame
            vb = obj.qrotate(obj.qconj(q),v);
            
            %Thruster forces + torques
            [Ft, Mt, vi0] = obj.thruster(thr,vb,w);
            
            %Slipstream velocities
            vstream = obj.slipstream(vi0);
            
            %Right Wing
            Frw = [0; 0; 0];
            Mrw = [0; 0; 0];
            for k = 1:length(obj.geom.lwing.span)
                xs = (obj.geom.rwing.pos{k} - obj.geom.cg)*1e-3; %position relative to CoM in meters
                ss = obj.geom.rwing.span(k)*1e-3; %span in meters
                cs = obj.geom.rwing.chord(k)*1e-3; %chord in meters
                As = ss*cs; %area in m^2
                cfs = obj.geom.rwing.flap(k)*1e-3; %flap chord in meters
                
                vs = vb + cross(w,xs) + [vstream.wing(k); 0; 0]; %velocity including slipstream
                as = atan2(vs(3),vs(1)); %angle of attack
                
                R = [cos(as) 0 -sin(as); 0 1 0; sin(as) 0 cos(as)]; %rotation from wind frame to body frame
                
                [Cl,Cd,Cm] = obj.flatPlateAirfoil(as,ail,cs,cfs,obj.geom.rwing.AR);
                
                Fsection = 0.5*obj.rho*As*(vs(1)*vs(1)+vs(3)*vs(3))*R*[-Cd; 0; -Cl];
                
                Frw = Frw + Fsection;
                Mrw = Mrw + cross(xs, Fsection) + 0.5*obj.rho*As*(vs(1)*vs(1)+vs(3)*vs(3))*cs*[0; Cm; 0];
            end
            
            %Left Wing
            Flw = [0; 0; 0];
            Mlw = [0; 0; 0];
            for k = 1:length(obj.geom.lwing.span)
                xs = (obj.geom.rwing.pos{k} - obj.geom.cg)*1e-3; %position relative to CoM in meters
                ss = obj.geom.rwing.span(k)*1e-3; %span in meters
                cs = obj.geom.rwing.chord(k)*1e-3; %chord in meters
                As = ss*cs; %area in m^2
                cfs = obj.geom.rwing.flap(k)*1e-3; %flap chord in meters
                
                vs = vb + cross(w,xs) + [vstream.wing(k); 0; 0]; %velocity including slipstream
                as = atan2(vs(3),vs(1)); %angle of attack
                
                R = [cos(as) 0 -sin(as); 0 1 0; sin(as) 0 cos(as)]; %rotation from wind frame to body frame
                
                [Cl,Cd,Cm] = obj.flatPlateAirfoil(as,-ail,cs,cfs,obj.geom.lwing.AR);
                
                Fsection = 0.5*obj.rho*As*(vs(1)*vs(1)+vs(3)*vs(3))*R*[-Cd; 0; -Cl];
                
                Flw = Flw + Fsection;
                Mlw = Mlw + cross(xs, Fsection) + 0.5*obj.rho*As*(vs(1)*vs(1)+vs(3)*vs(3))*cs*[0; Cm; 0];
            end
            
            %Horizontal Tail
            Fht = [0; 0; 0];
            Mht = [0; 0; 0];
            for k = 1:length(obj.geom.hTail.span)
                xs = (obj.geom.hTail.pos{k} - obj.geom.cg)*1e-3; %position relative to CoM in meters
                ss = obj.geom.hTail.span(k)*1e-3; %span in meters
                cs = obj.geom.hTail.chord(k)*1e-3; %chord in meters
                As = ss*cs; %area in m^2
                cfs = obj.geom.hTail.flap(k)*1e-3; %flap chord in meters
                
                vs = vb + cross(w,xs) + [vstream.hTail(k); 0; 0]; %velocity including slipstream
                as = atan2(vs(3),vs(1)); %angle of attack
                
                R = [cos(as) 0 -sin(as); 0 1 0; sin(as) 0 cos(as)]; %rotation from wind frame to body frame
                
                [Cl,Cd,Cm] = obj.flatPlateAirfoil(as,ele,cs,cfs,obj.geom.hTail.AR);
                
                Fsection = 0.5*obj.rho*As*(vs(1)*vs(1)+vs(3)*vs(3))*R*[-Cd; 0; -Cl];
                
                Fht = Fht + Fsection;
                Mht = Mht + cross(xs, Fsection) + 0.5*obj.rho*As*(vs(1)*vs(1)+vs(3)*vs(3))*cs*[0; Cm; 0];
            end
            
            %Vertical Tail
            Fvt = [0; 0; 0];
            Mvt = [0; 0; 0];
            for k = 1:length(obj.geom.vTail.span)
                xs = (obj.geom.vTail.pos{k} - obj.geom.cg)*1e-3; %position relative to CoM in meters
                ss = obj.geom.vTail.span(k)*1e-3; %span in meters
                cs = obj.geom.vTail.chord(k)*1e-3; %chord in meters
                As = ss*cs; %area in m^2
                cfs = obj.geom.vTail.flap(k)*1e-3; %flap chord in meters
                
                vs = vb + cross(w,xs) + [vstream.hTail(k); 0; 0]; %velocity including slipstream
                as = atan2(vs(2),vs(1)); %angle of attack
                
                R = [cos(as) -sin(as) 0; sin(as) cos(as) 0; 0 0 1]; %rotation from wind frame to body frame
                
                [Cl,Cd,Cm] = obj.flatPlateAirfoil(as,-rud,cs,cfs,obj.geom.vTail.AR);
                
                Fsection = 0.5*obj.rho*As*(vs(1)*vs(1)+vs(3)*vs(3))*R*[-Cd; -Cl; 0];
                
                Fvt = Fvt + Fsection;
                Mvt = Mvt + cross(xs, Fsection) + 0.5*obj.rho*As*(vs(1)*vs(1)+vs(3)*vs(3))*cs*[0; 0; Cm];
            end
            
            %Body
            Fb = [0; 0; 0];
            Mb = [0; 0; 0];
            for k = 1:length(obj.geom.body.span)
                xs = (obj.geom.body.pos{k} - obj.geom.cg)*1e-3; %position relative to CoM in meters
                ss = obj.geom.body.span(k)*1e-3; %span in meters
                cs = obj.geom.body.chord(k)*1e-3; %chord in meters
                As = ss*cs; %area in m^2
                cfs = obj.geom.vTail.flap(k)*1e-3; %flap chord in meters
                
                vs = vb + cross(w,xs) + [vstream.hTail(k); 0; 0]; %velocity including slipstream
                as = atan2(vs(2),vs(1)); %angle of attack
                
                R = [cos(as) -sin(as) 0; sin(as) cos(as) 0; 0 0 1]; %rotation from wind frame to body frame
                
                [Cl,Cd,Cm] = obj.flatPlateAirfoil(as,0,cs,cfs,obj.geom.vTail.AR);
                
                Fsection = 0.5*obj.rho*As*(vs(1)*vs(1)+vs(3)*vs(3))*R*[-Cd; -Cl; 0];
                
                Fb = Fb + Fsection;
                Mb = Mb + cross(xs, Fsection) + 0.5*obj.rho*As*(vs(1)*vs(1)+vs(3)*vs(3))*cs*[0; 0; Cm];
            end
            
            %Add everything together
            F = Ft + Frw + Flw + Fht + Fvt + Fb;
            M = Mt + Mrw + Mlw + Mht + Mvt + Mb;
        end
        
        function g = geomSetup(obj)
            g = struct();
            
            g.prop = struct();
            g.rwing = struct();
            g.lwing = struct();
            g.hTail = struct();
            g.vTail = struct();
            g.body = struct();
            
            g.cg = [-270 0 5.89]'; %position of CG relative to propeller
            
            g.prop.d = 254.0; %overall prop diameter (blade + hub)
            g.prop.dh = 12.0; %diameter of prop hub
            
            g.rwing.AR = 2.092;
            g.rwing.span = [57.26 58.72 59.29 62.58 62.58 63.25 61.90]';
            g.rwing.area = [14499.66 13977.71 13198.25 12974.40 11992.21 11122.51 9913.9]';
            g.rwing.chord = [253.4000  237.5800  222.6900  207.4200  191.7400  175.9700  160.2800]';
            g.rwing.flap = [0  102.0400   99.2300   96.3600   93.4100   90.4400   87.4900]';
            g.rwing.pos = {[-217.0700   33.3800         0]', ... 
                          [-211.3900   91.7600         0]', ...
                          [-215.4300  151.1400         0]', ...
                          [-212.1800  212.0100         0]', ...
                          [-208.8500  274.5500         0]', ...
                          [-205.5000  337.4200         0]', ...
                          [-202.2600  400.2100         0]'};
            
            g.lwing.AR = 2.092;
            g.lwing.span = [57.26 58.72 59.29 62.58 62.58 63.25 61.90]';
            g.lwing.area = [14499.66 13977.71 13198.25 12974.40 11992.21 11122.51 9913.9]';
            g.lwing.chord = [253.4000  237.5800  222.6900  207.4200  191.7400  175.9700  160.2800]';
            g.lwing.flap = [0  102.0400   99.2300   96.3600   93.4100   90.4400   87.4900]';
            g.lwing.pos = {[-217.0700   -33.3800         0]', ... 
                          [-211.3900   -91.7600         0]', ...
                          [-215.4300  -151.1400         0]', ...
                          [-212.1800  -212.0100         0]', ...
                          [-208.8500  -274.5500         0]', ...
                          [-205.5000  -337.4200         0]', ...
                          [-202.2600  -400.2100         0]'};
                      
            g.hTail.AR = 1.2421;
            g.hTail.span = [66.5500   51.0900   62.0200 66.5500   51.0900   62.0200]';
            g.hTail.area = [7241.76    7452.58    8892.67 7241.76    7452.58    8892.67]';
            g.hTail.chord = [117.4100  145.6400  143.1600 117.4100  145.6400  143.1600]';
            g.hTail.flap = [81.8700  111.4900  143.1600 81.8700  111.4900  143.1600]';
            g.hTail.pos = {[-711.2700   39.5300         0]', ...
                           [-719.7200   93.3300         0]', ...
                           [-715.1500  151.5700         0]', ...
                           [-711.2700   39.5300         0]', ...
                           [-719.7200   93.3300         0]', ...
                           [-715.1500  151.5700         0]'};
            
            g.vTail.AR = 1.2984;
            g.vTail.span = [46.3000   62.1900   69.3100   38.1000]';
            g.vTail.area = [9291.95  11997.38  11715.82   5298.95]';
            g.vTail.chord = [200.6200  192.9500  169.4900  141.1400]';
            g.vTail.flap = [120.3600  112.6900  103.4300  141.1400]';
            g.vTail.pos = {[-733.9400         0   23.5400]', ...
                           [-735.8100         0  -29.9600]', ...
                           [-744.1400         0  -94.6300]', ...
                           [-760.7600         0 -148.5600]'};
                       
            g.body.AR = 0.2279;
            g.body.span = [35.8700   35.8700   37.2600   37.2600]';
            g.body.area = [23018.14  23018.14  23910.11  23910.11]';
            g.body.chord = [641.7100  641.7100  641.7100  641.7100]';
            g.body.flap = [0     0     0     0]';
            g.body.pos = {[-205.0900         0   53.8100]', ...
                          [-205.0900         0   17.9400]', ...
                          [-205.0900         0  -18.6300]', ...
                          [-205.0900         0  -55.8900]'};        
        end
        
        function [Cl,Cd,Cm] = flatPlateAirfoil(obj,a,d,c,cf,AR)
            %Inputs:
            %a = angle of attack
            %d = flap deflection angle
            %c = chord
            %df = flap chord
            %AR = aspect ratio
            
            CL_a = 2*pi*AR/(AR+ 2*(AR+4)/(AR+2)); %Cl slope at zero AoA
            Kp = CL_a;
            Kv = pi;
            
            %Flap deflection calculations
            theta_f = acos(2*cf/c - 1);
            tau = 1 - (theta_f - sin(theta_f))/pi;
            eta = obj.eta(d);
            
            delta_CL = CL_a*tau*eta*d;
%             delta_CLmax = obj.ClmaxFactor(cf/c)*delta_CL;
            
            a0_eff = obj.fastNewton(@(x)obj.lift(Kp,Kv,delta_CL,-x), -delta_CL/Kp);
            
%             a_low = obj.lowAlphaLimit(AR);
%             a_high = obj.highAlphaLimit(AR);
%             
%             CLmax = Kp*sin(a_low)*(cos(a_low))^2 + Kv*sas(a_low)*cos(a_low) + delta_CLmax;
%             CLmin = Kp*sin(-a_low)*(cos(-a_low))^2 + Kv*sas(-a_low)*cos(-a_low) + delta_CLmax;
            
%             a_CLmax = obj.fastNewton(@(x)obj.lift(Kp,Kv,CLmax),a_low);
%             a_CLmin = obj.fastNewton(@(x)obj.lift(Kp,Kv,CLmin),-a_low);
            
            % High AoA parameters using equivalent flat plate method
%             c_eff = sqrt((c-cf)^2 + cf^2 + 2*(c-cf)*cf*cos(d));
%             gamma = asin(sin(d)*cf/c_eff);
            
            %Low-alpha model
            a_eff = a - a0_eff;
            CL1 = Kp*sin(a_eff)*(cos(a_eff))^2 + Kv*obj.sas(a_eff)*cos(a_eff);
            CD1 = Kp*cos(a_eff)*(sin(a_eff))^2 + Kv*obj.sas(a_eff)*sin(a_eff) + obj.Cd0;
            CM1 = -0.17*Kv*obj.sas(a_eff);
            
            %High-alpha model
%             CN2 = Cd90_eff*(
%             CT2 = 0.5*Cd0*cos(a_eff)
%             CL2 = 
%             CD2 = 
%             CM2 = 
            
            %Interpolation
            Cl = CL1;
            Cd = CD1;
            Cm = CM1;
            
        end
        
        function [F,M,Vi0] = thruster(obj,thr,vb,w)

            %thr command is prop angular velocity in rad/sec
            D = obj.geom.prop.d*1e-3; %prop diameter
            rprop = -obj.geom.cg*1e-3; %vector from CoM to prop
            
            vprop = vb + cross(w,rprop); %prop velocity
            vmag = sqrt(vprop'*vprop + 1e-8); %magnitude of velocity (padded to avoid badness near zero)
            va = sqrt(vprop(1)'*vprop(1) + 1e-8); %magnitude of axial velocity (padded to avoid badness near zero)
            vp = sqrt(vprop(2:3)'*vprop(2:3) + 1e-8); %magnitude of in-plane velocity (padded to avoid badness near zero)
            
            psi = atan2(vp,va); % azimuth angle of prop
            delta = atan2(vprop(3),vprop(2)); % in-plane angle of prop
            
            Jp = vmag*2*pi/(thr*D); %advance ratio
            
            %Evaluate forces and moments in propeller frame
            Cfx = CFx(Jp,psi);
            Cfy = CFy(Jp,psi);
            Cmx = CMx(Jp,psi);
            Cmy = CMy(Jp,psi);
            Cmz = CMz(Jp,psi);
            Fprop = obj.rho*(thr/(2*pi))^2*D^4*[Cfx; 0.5*Cfy; 0]; %no idea where these 1/2 factors come from...
            Mprop = obj.rho*(thr/(2*pi))^2*D^5*[Cmx; 0.5*Cmy; 0.5*Cmz];
            
            %Rotate into body frame
            R = [1     0          0;
                 0 -cos(delta)  sin(delta);
                 0 -sin(delta) -cos(delta)];
            F = R*Fprop;
            M = R*Mprop + cross(rprop,F);
            
            Vi0 = 0.5*0.2531*thr*D*tsqrt(Cfx);
            
            function s = tsqrt(x)
                t = obj.smoothstep(5*x);
                s = sqrt(t*x);
            end
            
            function c = CFx(x, y)
                %Polynomial fit to prop data from Khan and Nahon code
                %x = advance ratio
                %y = thruster azimuth angle
                
                p00 =      0.1533;
                p10 =     -0.1206;
                p01 =     0.06794;
                p20 =     -0.1636;
                p11 =     -0.1482;
                p02 =     -0.1958;
                p30 =     0.01551;
                p21 =     0.08877;
                p12 =      0.4111;
                p03 =      0.1819;
                p40 =     0.02015;
                p31 =    -0.04616;
                p22 =     0.04006;
                p13 =     -0.1711;
                p04 =    -0.05355;
                
                c = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 ...
                    + p30*x^3 + p21*x^2*y + p12*x*y^2 + p03*y^3 + p40*x^4 ...
                    + p31*x^3*y + p22*x^2*y^2 + p13*x*y^3 + p04*y^4;
            end
            
            function c = CFy(x, y)
                %Polynomial fit to prop data from Khan and Nahon code
                %x = advance ratio
                %y = thruster azimuth angle
                
                p00 =  -6.772e-06;
                p10 =   -0.003032;
                p01 =    0.001488;
                p20 =     0.00832;
                p11 =     0.02789;
                p02 =   -0.003636;
                p30 =   -0.006045;
                p21 =    0.003672;
                p12 =    -0.01236;
                p03 =    0.001903;
                
                c = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 ...
                    + p30*x^3 + p21*x^2*y + p12*x*y^2 + p03*y^3;
            end
            
            function c = CMx(x, y)
                %Polynomial fit to prop data from Khan and Nahon code
                %x = advance ratio
                %y = thruster azimuth angle
                
                p00 =   -0.009669;
                p10 =   -0.005389;
                p01 =    0.004559;
                p20 =      0.0207;
                p11 =    0.004111;
                p02 =   -0.007563;
                p30 =    0.007277;
                p21 =      -0.024;
                p12 =     0.00213;
                p03 =    0.003058;
                
                c = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 ...
                    + p30*x^3 + p21*x^2*y + p12*x*y^2 + p03*y^3;
            end
            
            function c = CMy(x, y)
                %Polynomial fit to prop data from Khan and Nahon code
                %x = advance ratio
                %y = thruster azimuth angle
                
                p00 =  -0.0005504;
                p10 =  -4.973e-05;
                p01 =     0.01133;
                p20 =    0.003695;
                p11 =     0.02263;
                p02 =    -0.03178;
                p30 =   -0.001043;
                p21 =    -0.03209;
                p12 =     0.03681;
                p03 =     0.02839;
                p40 =   -0.001924;
                p31 =    0.003947;
                p22 =     0.01706;
                p13 =    -0.02508;
                p04 =   -0.007967;
                
                c = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 ...
                    + p30*x^3 + p21*x^2*y + p12*x*y^2 + p03*y^3 + p40*x^4 ...
                    + p31*x^3*y + p22*x^2*y^2 + p13*x*y^3 + p04*y^4;
            end
            
            function c = CMz(x, y)
                %Polynomial fit to prop data from Khan and Nahon code
                %x = advance ratio
                %y = thruster azimuth angle
                
                p00 =   4.334e-05;
                p10 =   -0.000234;
                p01 =   -0.006584;
                p20 =      0.0081;
                p11 =     0.03073;
                p02 =     0.03908;
                p30 =    -0.02703;
                p21 =    -0.06918;
                p12 =    -0.03117;
                p03 =    -0.07128;
                p40 =     0.02886;
                p31 =     0.07995;
                p22 =    -0.03062;
                p13 =     0.06609;
                p04 =     0.05073;
                p50 =   -0.009777;
                p41 =    -0.03357;
                p32 =    0.009717;
                p23 =     0.01006;
                p14 =    -0.02981;
                p05 =    -0.01239;
                
                c = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 ...
                    + p30*x^3 + p21*x^2*y + p12*x*y^2 + p03*y^3 + p40*x^4 ...
                    + p31*x^3*y + p22*x^2*y^2 + p13*x*y^3 + p04*y^4 ...
                    + p50*x^5 + p41*x^4*y + p32*x^3*y^2 + p23*x^2*y^3 ...
                    + p14*x*y^4 + p05*y^5;
            end
            
        end
        
        function vstream = slipstream(obj, vi0)
            %Implements the propeller slipstream model from Khan and Nahon (2015)
            
            rp = obj.geom.prop.d/(2e3);
            rh = obj.geom.prop.dh/(2e3);
            
            x0 = 1.528*rp; %eq. 11
            R0 = 0.74*rp;  %eq. 12
            D0 = 2*R0;
            rmax0 = 0.67*(rp-rh);
            v0 = 2*vi0*(1.46/1.59); %efflux velocity
            
            
            vstream.wing = zeros(length(obj.geom.rwing.pos),1);
            for k = 1:length(obj.geom.rwing.pos)
                vstream.wing(k) = vfit(obj.geom.rwing.pos{k});
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
        
        function [s, ds] = smoothstep(obj,x)
            %C2 interpolation from 0 to 1
            if x < 0
                s = 0;
                ds = 0;
            elseif x < 1
                s = x*x*x*(x*(x*6 - 15) + 10); %6*x^5 - 15*x^4 + 10*x^3;
                ds = x*x*(x*(x*30 - 60) + 30); %30*x^4 - 60*x^3 + 30*x^2;
            else
                s = 1;
                ds = 0;
            end
        end
        
        function [s, ds] = sas(obj,x)
            %A smooth approximation of the function sin(x)*abs(sin(x))
            
            a = -0.105;
            b =  5.682;
            w1 = 3.944;
            w2 = 0.1127;
            
            s = a*sin(w1*x) + b*sin(w2*x);
            ds = a*w1*cos(w1*x) + b*w2*cos(w2*x);
        end
        
        function f = ClmaxFactor(obj, r)
            %Cubic polynomial fit to CLmax reduction factor as a function
            %of flap chord ratio (see eq. 3.16 in Khan)
            
            p1 =       1.343;
            p2 =      -1.916;
            p3 =     -0.4182;
            p4 =      0.9971;
            
            f = r*(p3 + r*(p2 + r*p1)) + p4;
            
        end
        
        function [e, de] = eta(obj, x)
            %Smooth curve fit to viscosity effect factor in flap deflection
            %lift calculation (eq. 3.14 in Khan)
            
            a =       0.081;
            b =     -0.3168;
            c =      0.6025;
            
            e = a*x.*x + b*sqrt(x.*x + .01) + c;
            de = 2*a*x - b*x./(sqrt(x.*x + .01));
        end
        
        function [Cl, dCl] = lift(obj,Kp,Kv,dCl,a)
            s = sin(a);
            c = cos(a);
            [sa, dsa] = obj.sas(a);
            
            %CL1 = Kp*sin(a_eff)*(cos(a_eff))^2 + Kv*obj.sas(a_eff)*cos(a_eff);
            
            Cl = Kp*s*c*c + Kv*sa*c - dCl;
            dCl = Kp*(c*c*c - 2*s*s*c) + Kv*(dsa*c - sa*s);
        end
        
        function a_low = lowAlphaLimit(obj,AR)
            %Valid between .5 and 6
            
            p1 =   -0.005595;
            p2 =     0.07587;
            p3 =     -0.3509;
            p4 =      0.7223;
            
            a_low = p4 + AR*(p3 + AR*(p2 + AR*p1));
        end
        
        function a_high = highAlphaLimit(obj,AR)
            %Valid between .5 and 6
            
            p1 =      -0.006;
            p2 =     0.08564;
            p3 =     -0.4097;
            p4 =       1.025;
       
            a_high = p4 + AR*(p3 + AR*(p2 + AR*p1));
        end
        
        function x = fastNewton(obj, f, x0)
            tol = 1e-7;
            x = x0;
            [e, de] = f(x0);
            while abs(e) > tol
                x = x + de\e;
                [e, de] = f(x);
            end
        end
    end
    
end

