function [CN,CL,CD,CM,Reg] = flappedAirfoil_LAR_v1_FlatPlate(i,alp,cf,c,def,AR,Cd0,Cd90,alp0)
k_CN = 1;

%----------------- STALL ANGLES AND HIGH AOA ANGLES -----------------%
% AR_pre = [0.5,0.75,1,1.25,1.5,1.75,2, 3,4,6];
% alpStall_pre = [44,37,34,22,20,18,18, 18,17,17]*pi/180;
% alpLim_pre = [35,33,28,20,15,14,13, 13,13,13]*pi/180; % By Torres    
AR_pre = [0.1666,0.333,0.4,0.5,1,1.25, 2, 3, 4, 6];
%alpLim_pre = [32,36,34,38,32,21,16,11,10,8]*pi/180; % Alp_Limit
LowAlpEnd_pre = [32,36,34,38,32,21,16,11,10,8]*pi/180; % Alp_Limit
HighAlpStart_pre = [40,60,55,56,40, 29, 28, 24, 22, 20]*pi/180;

LowAlpEnd = interp1(AR_pre,LowAlpEnd_pre,AR);
HighAlpStart = interp1(AR_pre, HighAlpStart_pre, AR);

%----------------------- LIFT CURVE SLOPE ---------------------------%
CLAlp_Method = 'McCormick';

switch CLAlp_Method
    case 'McCormick' % From McCormick
        ClAlp = 2*pi;
        CLAlp = ClAlp*AR/(AR + 2*(AR+4)/(AR+2));

    case 'Anderson' % FROM ANDERSON (THROUGH TORRES AND MUELLER)
        ClAlp = 5.3743;
        tau = 0.3; %%%%
        CLAlp = ClAlp/(1 + (1 + tau)*ClAlp/(pi*AR));

    case 'Polhamus' % FROM LOWRY & POLHAMUS (THROUGH TORRES)
        ClAlp = 2*pi;
        eta = ClAlp/(2*pi);
        CLAlp = 2*pi*AR/(2 + sqrt((AR/eta)^2*1 + 4));

    case 'Hoerner' % FROM HOERNER
        CLAlp = 1/(36.5/AR + 2*AR)*180/pi;
end

% -------- POTENTIAL LIFT AND VORTEX LIFT COEFFICIENTS --------------%
Kp = CLAlp;
Kv = pi;

%--------------------- FLAP DEFLECTION ------------------------------%
if cf/c == 1
    alp0_eff = alp0;

    LowAlpEnd_P = LowAlpEnd;
    LowAlpEnd_N = -LowAlpEnd;
    LowAlpStart_P = LowAlpEnd_N + pi;
    LowAlpStart_N = LowAlpEnd_P - pi;

    HighAlpStart_P = HighAlpStart;
    HighAlpEnd_P = pi - HighAlpStart_P;
    HighAlpStart_N = -HighAlpStart_P;
    HighAlpEnd_N = -HighAlpEnd_P;

    alp = alp + def;
    c_eff = c;
    gamma = 0;                           

else
    the_f = acos(2*cf/c - 1);
    tau = 1 - (the_f - sin(the_f))/pi;

    cf2c_pre = [0 0.2 0.4 0.6 0.8 1];
    delCLmax2delCL_pre = [1 0.83 0.65 0.3 0.15 0];
    delCLmax2delCL = interp1(cf2c_pre,delCLmax2delCL_pre,cf/c);

    beta_pre = [0 10 20 30 40 50 60 70]*pi/180;    
    %eta_pre = [0.81*0.7 0.8*0.7 0.685*0.7 0.535 0.455 0.405 0.37 0.345];
    eta_pre = [0.582, 0.536, 0.494, 0.456, 0.42, 0.388, 0.359, 0.333];
    eta = interp1(beta_pre, eta_pre, abs(def));

    delCL = CLAlp*tau*eta*def;
    delCLmax = delCLmax2delCL*delCL;

    % For Linear CL vs Alp curve
    % alp0_eff1 = alp0 - delCL/CLAlp;                
    % CLmaxP1 = CLAlp*(alpStallP - alp0) + delCLmax;
    % CLmaxN1 = CLAlp*(alpStallN - alp0) + delCLmax;
    % alpStallP_eff1 = alp0_eff1 + CLmaxP1/CLAlp;
    % alpStallN_eff1 = alp0_eff1 + CLmaxN1/CLAlp;

    % For Non-Linear CL vs Alp curve
    %fun = @(x) Kp*sin(0-x)*(cos(0-x))^2 + Kv*abs(sin(0-x))*sin(0-x)*cos(0-x) - delCL;
    x0 = alp0;              
    alp0_eff = fzero(@myfun_alp0_eff,x0,[],Kp,Kv,delCL);

    a_eff = LowAlpEnd - alp0;        
    CLmaxP = Kp*sin(a_eff)*(cos(a_eff))^2 + Kv*abs(sin(a_eff))*sin(a_eff)*cos(a_eff) + delCLmax;
    a_eff = -LowAlpEnd - alp0;
    CLmaxN = Kp*sin(a_eff)*(cos(a_eff))^2 + Kv*abs(sin(a_eff))*sin(a_eff)*cos(a_eff) + delCLmax;

    %fun = @(x) Kp*sin(x-alp0_eff)*(cos(x-alp0_eff))^2 + Kv*abs(sin(x-alp0_eff))*sin(x-alp0_eff)*cos(x-alp0_eff) - CLmaxP;
    x0 = LowAlpEnd;
    LowAlpEnd_P = fzero(@myfun_alpStall_eff,x0,[],Kp,Kv,alp0_eff,CLmaxP);
    %alpStallP_eff = fzero(@myfun_alpStall_eff,x0,[],Kp,Kv,alp0_eff,CLmaxP);

    %fun = @(x) Kp*sin(x-alp0_eff)*(cos(x-alp0_eff))^2 + Kv*abs(sin(x-alp0_eff))*sin(x-alp0_eff)*cos(x-alp0_eff) - CLmaxN;
    x0 = -LowAlpEnd;
    LowAlpEnd_N = fzero(@myfun_alpStall_eff,x0,[],Kp,Kv,alp0_eff,CLmaxN);
    %alpStallN_eff = fzero(@myfun_alpStall_eff,x0,[],Kp,Kv,alp0_eff,CLmaxN);

    LowAlpStart_P = LowAlpEnd_N + pi;
    LowAlpStart_N = LowAlpEnd_P - pi;

    % High AoA parameters using equivalent flat plate method
    c_eff = sqrt((c-cf)^2 + cf^2 + 2*(c-cf)*cf*cos(abs(def)));
    gamma = asin(sin(def)*cf/c_eff);

    HighAlpStart_P = HighAlpStart;
    HighAlpEnd_P = pi - HighAlpStart_P;
    HighAlpStart_N = -HighAlpStart_P;
    HighAlpEnd_N = -HighAlpEnd_P;
end

    %% Calculate Aerodynamic coefficients based on regime    
    if alp >= -pi && alp <= LowAlpStart_N
        Reg = -4;
        a_eff = alp - alp0_eff + pi;
        
        CL = Kp*sin(a_eff)*(cos(a_eff))^2 + Kv*abs(sin(a_eff))*sin(a_eff)*cos(a_eff);        
        CD = Cd0 + abs(Kp*abs(sin(a_eff))*sin(a_eff)*cos(a_eff) + Kv*(sin(a_eff))^3);        
        CN = -(Kp*sin(a_eff)*abs(cos(a_eff)) + Kv*abs(sin(a_eff))*sin(a_eff));        
        %CM = (0.42-0.25)*Kv*abs(sin(a_eff))*sin(a_eff);
        CM = -(0.25-0.175*(1-abs(a_eff-pi)/(pi/2)))*CN;
    
    % ------------------------------------------------- %    
    elseif alp > LowAlpStart_N && alp < HighAlpEnd_N
        Reg = -3;
        % Linear End
        a_eff = LowAlpStart_N - alp0_eff + pi;
        
        CL1 = Kp*sin(a_eff)*(cos(a_eff))^2 + Kv*abs(sin(a_eff))*sin(a_eff)*cos(a_eff);
        CD1 = Cd0 + abs(Kp*abs(sin(a_eff))*sin(a_eff)*cos(a_eff) + Kv*(sin(a_eff))^3);
        CN1 = -(Kp*sin(a_eff)*cos(a_eff) + Kv*abs(sin(a_eff))*sin(a_eff));
        %CM1 = (0.42-0.25)*Kv*abs(sin(a_eff))*sin(a_eff);        
        %CD1 = Cd0 + abs(CL1)*abs(tan(a_eff));
        %CN1 = -sqrt(CL1^2 + CD1^2);        
        CM1 = -(0.25-0.175*(1-abs(a_eff-pi)/(pi/2)))*CN1;           
        
        % NonLinear End
        a_eff = HighAlpEnd_N + gamma; %- alp0_eff;
        
        % Determining effective Cd90
        if cf/c == 1 || cf/c == 0
            Cd90_eff = Cd90;
        elseif sign(HighAlpEnd_N) == 1
            Cd90_eff = -1.29630E-05*(def*180/pi)^2 + 3.61111E-03*(def*180/pi) + Cd90;
        else
            Cd90_eff = -1.29630E-05*(def*180/pi)^2 - 3.61111E-03*(def*180/pi) + Cd90;
        end
        
        % Hoerner Model with Lindenberg Correction
        CN2 = -k_CN*Cd90_eff*(1/(0.56 + 0.44*sin(abs(a_eff))) - 0.41*(1-exp(-17/AR)))*sin(abs(a_eff));
        CT2 = 0.5*Cd0*cos(a_eff);
        CL2 = CN2*cos(a_eff) - CT2*sin(a_eff);
        CD2 = CN2*sin(a_eff) + CT2*cos(a_eff);
        CM2 = -(0.25-0.175*(1-abs(a_eff)/(pi/2)))*CN2;
        
        % Interpolation
        CN = interp1([LowAlpStart_N,HighAlpEnd_N],[CN1,CN2],alp);
        CL = interp1([LowAlpStart_N,HighAlpEnd_N],[CL1,CL2],alp);
        CD = interp1([LowAlpStart_N,HighAlpEnd_N],[CD1,CD2],alp);
        CM = interp1([LowAlpStart_N,HighAlpEnd_N],[CM1,CM2],alp);
        
    % ------------------------------------------------- %    
    elseif alp >= HighAlpEnd_N && alp <= HighAlpStart_N
        Reg = -2;
        % Equivalent flat plate
        a_eff = alp + gamma; %- alp0_eff;
        
        % Determining effective Cd90
        if cf/c == 1 || cf/c == 0
            Cd90_eff = Cd90;
        elseif sign(alp) == 1
            Cd90_eff = -1.29630E-05*(def*180/pi)^2 + 3.61111E-03*(def*180/pi) + Cd90;
        else
            Cd90_eff = -1.29630E-05*(def*180/pi)^2 - 3.61111E-03*(def*180/pi) + Cd90;
        end
        
        % Hoerner Model
        CN = -k_CN*Cd90_eff*(1/(0.56 + 0.44*sin(abs(a_eff))) - 0.41*(1-exp(-17/AR)))*sin(abs(a_eff));
        CT = 0.5*Cd0*cos(a_eff);
        CL = CN*cos(a_eff) - CT*sin(a_eff);
        CD = CN*sin(a_eff) + CT*cos(a_eff);
        CM = -(0.25-0.175*(1-abs(a_eff)/(pi/2)))*CN; 
        
    % ------------------------------------------------- %
    elseif alp > HighAlpStart_N && alp < LowAlpEnd_N
        Reg = -1;
        % Linear End
        a_eff = LowAlpEnd_N - alp0_eff;
        CL1 = Kp*sin(a_eff)*(cos(a_eff))^2 + Kv*abs(sin(a_eff))*sin(a_eff)*cos(a_eff);
        CD1 = Cd0 + abs(Kp*abs(sin(a_eff))*sin(a_eff)*cos(a_eff) + Kv*(sin(a_eff))^3);
        CN1 = Kp*sin(a_eff)*cos(a_eff) + Kv*abs(sin(a_eff))*sin(a_eff);
        CM1 = -(0.42-0.25)*Kv*abs(sin(a_eff))*sin(a_eff);          
        
        % NonLinear End
        a_eff = HighAlpStart_N + gamma; %- alp0_eff;

        if cf/c == 1 || cf/c == 0
            Cd90_eff = Cd90;
        elseif sign(HighAlpStart_N) == 1
            Cd90_eff = -1.29630E-05*(def*180/pi)^2 + 3.61111E-03*(def*180/pi) + Cd90;
        else
            Cd90_eff = -1.29630E-05*(def*180/pi)^2 - 3.61111E-03*(def*180/pi) + Cd90;
        end

        CN2 = -k_CN*Cd90_eff*(1/(0.56 + 0.44*sin(abs(a_eff))) - 0.41*(1-exp(-17/AR)))*sin(abs(a_eff));
        CT2 = 0.5*Cd0*cos(a_eff);    
        CL2 = CN2*cos(a_eff) - CT2*sin(a_eff);
        CD2 = CN2*sin(a_eff) + CT2*cos(a_eff);
        CM2 = -(0.25-0.175*(1-abs(a_eff)/(pi/2)))*CN2;        
        
        % Interpolation
        CN = interp1([LowAlpEnd_N,HighAlpStart_N],[CN1,CN2],alp);
        CL = interp1([LowAlpEnd_N,HighAlpStart_N],[CL1,CL2],alp);
        CD = interp1([LowAlpEnd_N,HighAlpStart_N],[CD1,CD2],alp);
        CM = interp1([LowAlpEnd_N,HighAlpStart_N],[CM1,CM2],alp); 
        
    % ------------------------------------------------- %        
    elseif alp >= LowAlpEnd_N && alp <= LowAlpEnd_P
        Reg = 0;
        a_eff = alp - alp0_eff;
        CL = Kp*sin(a_eff)*(cos(a_eff))^2 + Kv*abs(sin(a_eff))*sin(a_eff)*cos(a_eff);        
        CD = Cd0 + abs(Kp*abs(sin(a_eff))*sin(a_eff)*cos(a_eff)) + Kv*abs((sin(a_eff))^3);
        CN = Kp*sin(a_eff)*cos(a_eff) + Kv*abs(sin(a_eff))*sin(a_eff);
        CM = -(0.42-0.25)*Kv*abs(sin(a_eff))*sin(a_eff);                    
        
    % ------------------------------------------------- %              
    elseif alp > LowAlpEnd_P && alp < HighAlpStart_P
        Reg = 1;
        % Linear End
        a_eff = LowAlpEnd_P - alp0_eff;
        CL1 = Kp*sin(a_eff)*(cos(a_eff))^2 + Kv*abs(sin(a_eff))*sin(a_eff)*cos(a_eff);
        CD1 = Cd0 + abs(Kp*abs(sin(a_eff))*sin(a_eff)*cos(a_eff) + Kv*(sin(a_eff))^3);
        CN1 = Kp*sin(a_eff)*cos(a_eff) + Kv*abs(sin(a_eff))*sin(a_eff);
        CM1 = -(0.42-0.25)*Kv*abs(sin(a_eff))*sin(a_eff);    
                
        % NonLinear End
        a_eff = HighAlpStart_P + gamma; %- alp0_eff;

        if cf/c == 1 || cf/c == 0
            Cd90_eff = Cd90;
        elseif sign(HighAlpStart_P) == 1
            Cd90_eff = -1.29630E-05*(def*180/pi)^2 + 3.61111E-03*(def*180/pi) + Cd90;
        else
            Cd90_eff = -1.29630E-05*(def*180/pi)^2 - 3.61111E-03*(def*180/pi) + Cd90;
        end

        CN2 = k_CN*Cd90_eff*(1/(0.56 + 0.44*sin(abs(a_eff))) - 0.41*(1-exp(-17/AR)))*sin(abs(a_eff));
        CT2 = 0.5*Cd0*cos(a_eff);
        CL2 = CN2*cos(a_eff) - CT2*sin(a_eff);
        CD2 = CN2*sin(a_eff) + CT2*cos(a_eff);
        CM2 = -(0.25-0.175*(1-abs(a_eff)/(pi/2)))*CN2;                       
        
        % Interpolation
        CN = interp1([LowAlpEnd_P,HighAlpStart_P],[CN1,CN2],alp);
        CL = interp1([LowAlpEnd_P,HighAlpStart_P],[CL1,CL2],alp);
        CD = interp1([LowAlpEnd_P,HighAlpStart_P],[CD1,CD2],alp);
        CM = interp1([LowAlpEnd_P,HighAlpStart_P],[CM1,CM2],alp);
        
    % ------------------------------------------------- %
    elseif alp >= HighAlpStart_P && alp <= HighAlpEnd_P
        Reg = 2;
        % Equivalent flat plate
        a_eff = alp + gamma; %- alp0_eff;
        
        % Determining effective Cd90
        if cf/c == 1 || cf/c == 0
            Cd90_eff = Cd90;
        elseif sign(alp) == 1
            Cd90_eff = -1.29630E-05*(def*180/pi)^2 + 3.61111E-03*(def*180/pi) + Cd90;
        else
            Cd90_eff = -1.29630E-05*(def*180/pi)^2 - 3.61111E-03*(def*180/pi) + Cd90;
        end
        
        % Hoerner Model
        CN = k_CN*Cd90_eff*(1/(0.56 + 0.44*sin(abs(a_eff))) - 0.41*(1-exp(-17/AR)))*sin(abs(a_eff));
        CT = 0.5*Cd0*cos(a_eff);
        CL = CN*cos(a_eff) - CT*sin(a_eff);
        CD = CN*sin(a_eff) + CT*cos(a_eff);
        CM = -(0.25-0.175*(1-abs(a_eff)/(pi/2)))*CN;
        
    % ------------------------------------------------- %    
    elseif alp > HighAlpEnd_P && alp < LowAlpStart_P
        Reg = -3;
        % Linear End
        a_eff = LowAlpStart_P - alp0_eff - pi;
        CL1 = Kp*sin(a_eff)*(cos(a_eff))^2 + Kv*abs(sin(a_eff))*sin(a_eff)*cos(a_eff);
        CD1 = Cd0 + abs(Kp*abs(sin(a_eff))*sin(a_eff)*cos(a_eff) + Kv*(sin(a_eff))^3);
        CN1 = -(Kp*sin(a_eff)*cos(a_eff) + Kv*abs(sin(a_eff))*sin(a_eff));
        CM1 = -(0.25-0.175*(1-abs(a_eff-pi)/(pi/2)))*CN1;
        %CM1 = (0.42-0.25)*Kv*abs(sin(a_eff))*sin(a_eff);          
        
        % NonLinear End
        a_eff = HighAlpEnd_P + gamma; %- alp0_eff;
        
        % Determining effective Cd90
        if cf/c == 1 || cf/c == 0
            Cd90_eff = Cd90;
        elseif sign(HighAlpEnd_P) == 1
            Cd90_eff = -1.29630E-05*(def*180/pi)^2 + 3.61111E-03*(def*180/pi) + Cd90;
        else
            Cd90_eff = -1.29630E-05*(def*180/pi)^2 - 3.61111E-03*(def*180/pi) + Cd90;
        end
        
        % Hoerner Model        
        CN2 = k_CN*Cd90_eff*(1/(0.56 + 0.44*sin(abs(a_eff))) - 0.41*(1-exp(-17/AR)))*sin(abs(a_eff));
        CT = 0.5*Cd0*cos(a_eff);
        CL2 = CN2*cos(a_eff) - CT*sin(a_eff);
        CD2 = CN2*sin(a_eff) + CT*cos(a_eff);
        CM2 = -(0.25-0.175*(1-abs(a_eff)/(pi/2)))*CN2;      
        
        % Interpolation
        CN = interp1([LowAlpStart_P,HighAlpEnd_P],[CN1,CN2],alp);
        CL = interp1([LowAlpStart_P,HighAlpEnd_P],[CL1,CL2],alp);
        CD = interp1([LowAlpStart_P,HighAlpEnd_P],[CD1,CD2],alp);
        CM = interp1([LowAlpStart_P,HighAlpEnd_P],[CM1,CM2],alp); 

    % ------------------------------------------------- %
    else
        Reg = 4;
        a_eff = alp - alp0_eff - pi;
        CL = Kp*sin(a_eff)*(cos(a_eff))^2 + Kv*abs(sin(a_eff))*sin(a_eff)*cos(a_eff);        
        CD = Cd0 + abs(Kp*abs(sin(a_eff))*sin(a_eff)*cos(a_eff) + Kv*(sin(a_eff))^3);
        CN = -(Kp*sin(a_eff)*cos(a_eff) + Kv*abs(sin(a_eff))*sin(a_eff));
        CM = -(0.25-0.175*(1-abs(a_eff-pi)/(pi/2)))*CN;
        %CM = (0.42-0.25)*Kv*abs(sin(a_eff))*sin(a_eff);    
    end     