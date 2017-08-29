function [f, df] = inverse_pzt_model_fun(obj, t, x, u)

xin = [t; x; u];
[f,df] = inverse_pzt_model(obj,xin);

df_fd = zeros(size(df));
step = 1e-6; %sqrt(eps(max(abs(xin))))
dxin = step*eye(length(xin));
for k = 1:length(xin)
    df_fd(:,k) = (inverse_pzt_model(obj,xin+dxin(:,k)) - ...
        inverse_pzt_model(obj,xin-dxin(:,k)))/(2*step);
end


disp('Inverse pzt model error:');
disp(max(abs(df_fd(:)-df(:))));
% max(abs(df(:)))


end

function [V, dV] = inverse_pzt_model(obj, xin)

nq = obj.getNumPositions();
nu = obj.getNumInputs();

t = xin(1);
x = xin(2:2*nq+1);
u = xin(2*nq+2:2*nq+1+nu);

act_joints = obj.getActuatedJoints();
qact = x(act_joints);     % convert back to SI
orientation = [1; 1; -1; -1; 1; 1; -1; -1];

% Drive Parameters
Vb = 200;
Vg = 0;

% Geometry
tact = 325e-3;          %m (actuator thickness)
tpzt = 135e-3;          %m (pzt thickness)
wn = 3.5;            %m (nominal width)
wr = 1.5;               %m (width ratio)
lext = 1;            %m (extension length)
lact = 9.2;          %m (actuator length)

% PZT internal Stress as a function of field (see N.T. Jafferis et al. 2015)
f31_eff_b = 24e-3;         %Pa*m/V (in blocked case)
f31_eff_f = 32e-3;         %Pa*m/V (in free case)

% PZT Young's modulus as a function of strain (see N.T. Jafferis et al. 2015)
Emin = 38.5e3;          %sigmoid min (Pa)
Emax = 81e3;            %sigmoid max (Pa)
Eave = 40e3;            % average modulus (Pa)
A0 = -0.3496;           %constant on cosine fit for Eave (1/m)
A1 = 0.3291;            %multiplier on cosine fit for Eave (1/m)
omg = 2.047;             %freq on cosine fit for Eave (unclear units)

% Material properties (see N.T. Jafferis et al. 2015)
Ecf = 340e3;             % Pa (CF young's modulus)

%  ----------- fillin in derived actuator properties -------------
tcf = tact - 2*tpzt;            % cf thickness
lr = lext/lact;                 % length ratio
eps = (Vb -Vg)/tpzt;        % total field

% geometric factor
GF_num = 8*(1-wr)^3*(1 + 2*lr);
GF_den = -6*(wr - 1)*(-3 + 4*lr*(wr - 1)+2*wr) + 3*(-2 + 2*lr*(wr-1)+wr)^2*log((2-wr)/wr);
GF = GF_num/GF_den;

% free deflection
df_den = (1/3)*Eave*tpzt*(1.5*tcf^2 + 3*tcf*tpzt + 2*tpzt^2) + Ecf*tcf^3/12;
df_num = f31_eff_f*tpzt*lact^2*eps*(tpzt+tcf);
df = 0.25*(1 + 2*lr)*(df_num/df_den);

E = Emin - (Emax - Emin)*(A0 + A1*cos(omg*qact));

kact = 0.85*3*(wn/lact^3)*GF*(((1/3)*E*tpzt*(1.5*tcf^2 + ...
    3*tcf*tpzt + 2*tpzt^2)+ Ecf*tcf^3/12)/(1 + 2*lr));

Ct = (0.75/lact)*(f31_eff_b + (qact/df)*(f31_eff_f - f31_eff_b))*wn*(tpzt+tcf)*GF;
Cb = (0.75/lact)*(f31_eff_b + (-qact/df)*(f31_eff_f - f31_eff_b))*wn*(tpzt+tcf)*GF;

V = (diag(Ct) + diag(Cb))\(Ct*Vb + Cb*Vg - diag(orientation)*(u + diag(kact)*qact));

% Derivative Stuff

% d/dx
dVdx = zeros(numel(u), numel(x));

dCt_dq = (0.75/(lact*df))*(f31_eff_f - f31_eff_b)*wn*(tpzt+tcf)*GF;
dCb_dq = -(0.75/(lact*df))*(f31_eff_f - f31_eff_b)*wn*(tpzt+tcf)*GF;

% dCt_dq + dCb_dq

dE_dq = (Emax - Emin)*(A1*omg*sin(omg*qact));
dkact_dq = dE_dq*(0.85*(wn/lact^3)*GF*tpzt*(1.5*tcf^2 + 3*tcf*tpzt + 2*tpzt^2))/(1+2*lr); 


for i = 1:numel(qact)
    d1dq = Vb*dCt_dq/(Ct(i) + Cb(i));
    d2dq = Vg*dCb_dq/(Ct(i) + Cb(i));
    d4dq = -orientation(i)*(kact(i) + dkact_dq(i)*qact(i))/(Ct(i) + Cb(i));
    dVdx(i, act_joints(i)) = d1dq + d2dq + d4dq;
end

%d/du
dVdu = -(diag(Ct) + diag(Cb))\diag(orientation);
dV = [zeros(numel(V), 1), dVdx, dVdu];
end