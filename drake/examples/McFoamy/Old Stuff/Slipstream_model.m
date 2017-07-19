function [VProp_Axial, VProp_Swirl] = Slipstream_model(RPM, Vi0_avg, Prop_Geom, Geometry)
%Rh = 0.006;
%Rp = 0.127;
Rh = Prop_Geom(1)/2*1e-3;
Rp = Prop_Geom(2)/2*1e-3;

DesX = abs(Geometry(5,5:end))*1e-3;
DesR = sqrt(Geometry(6,5:end).^2 + Geometry(7,5:end).^2)*1e-3;

VProp_Axial = zeros(1,length(DesX));
VProp_Swirl = zeros(1,length(DesX));

% Equations used from propeller slipstream paper
x0 = 1.528*Rp;   % Eq. 11
R0 = 0.74*Rp;   % Eq. 12
D0 = 2*R0;
Rm0 = 0.67*(R0 - Rh);
V0 = 2*Vi0_avg*(1.46/1.59); % Efflux velocity from Vi0_avg

for i = 1:length(DesX)
    if DesX(i) < x0
        VProp_Axial(i) = Vi0_avg*(1 + (DesX(i)/Rp)/sqrt(1 + (DesX(i)/Rp)^2));  % Eq. 1        

    elseif DesX(i) < 1.7*D0
        Vmax = V0*(1.24 - 0.0765*(DesX(i)-x0)/D0);  % Eq. 13
        Rm = Rm0*(1 - 0.1294*(DesX(i) - x0)/D0);    % Eq. 14
        VProp_Axial(i) = Vmax*exp(-((DesR(i)-Rm)/(0.8839*Rm0 + 0.1326*(DesX(i)-x0-R0)))^2);   % Eq. 15
    
    elseif DesX(i) < 4.25*D0
        Vmax = V0*(1.37 - 0.1529*(DesX(i)-x0)/D0);  % Eq. 16
        Rm = Rm0*(1.3 - 0.3059*(DesX(i) - x0)/D0);  % Eq. 17
        VProp_Axial(i) = Vmax*exp(-((DesR(i)-Rm)/(0.5176*Rm0 + 0.2295*(DesX(i)-x0-R0)))^2);   % Eq. 18
    
    else
        Vmax = V0*(0.89 - 0.04*(DesX(i)-x0)/D0);    % Eq. 19
        VProp_Axial(i) = Vmax*exp(-(DesR(i)/(0.2411*(DesX(i)-x0)))^2);    % Eq. 21
        
    end
    
    if VProp_Axial(i) < 0.01
        VProp_Axial(i) = 0; 
    end
end

% Calculate the radial component on Wing only
% DesR_Wing = DesR(1:6);
% VProp_Rad = zeros(size(DesR_Wing));
%
% r = [0,0.35759,0.41593,...
%     0.47616,0.5345,0.59849,0.65684,0.71518,0.77353,0.83375,0.89021]*Rp;
% 
% Vtip = Rp*RPM*2*pi/60;
% Vt = [0,0.035691875,...
%     0.031344793,0.02951442,0.03203119,0.031344793,0.034319141,0.032946346,0.029285662,0.021964231,0.017845968]*Vtip;
% 
% for i = 1:length(DesR_Wing)
%     if DesR_Wing(i) <= 0.89*Rp
%         VProp_Rad(i) = interp1(r,Vt,DesR_Wing(i));
%     else
%         VProp_Rad(i) = 0;
%     end
% end

