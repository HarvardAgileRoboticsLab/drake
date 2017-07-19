function [AilDefOut,ElevDefOut,RudDefOut,ThrRPMOut] = Command_saturation(AilDef_max,ElevDef_max,RudDef_max,ThrRPM_max,Controls)

AilDefIn = Controls(1);
ElevDefIn = Controls(2);
RudDefIn = Controls(3);
ThrRPMIn = Controls(4);


if AilDefIn >= -AilDef_max && AilDefIn <= AilDef_max
    AilDefOut = AilDefIn;
elseif AilDefIn < -AilDef_max
    AilDefOut = -AilDef_max;
else
    AilDefOut = AilDef_max;
end

if ElevDefIn >= -ElevDef_max && ElevDefIn <= ElevDef_max
    ElevDefOut = ElevDefIn;
elseif ElevDefIn < -ElevDef_max
    ElevDefOut = -ElevDef_max;
else
    ElevDefOut = ElevDef_max;
end

if RudDefIn >= -RudDef_max && RudDefIn <= RudDef_max
    RudDefOut = RudDefIn;
elseif RudDefIn < -RudDef_max
    RudDefOut = -RudDef_max;
else
    RudDefOut = RudDef_max;
end

if ThrRPMIn >= 0 && ThrRPMIn <= ThrRPM_max
    ThrRPMOut = ThrRPMIn;
elseif ThrRPMIn < 0
    ThrRPMOut = 0;
else
    ThrRPMOut = ThrRPM_max;
end