function [f, df] = ipzt_fun(obj, act, t, x, u)

xin = [t; x; u];
[f,df] = ipzt(obj,act, xin);

df_fd = zeros(size(df));
step = 1e-6; %sqrt(eps(max(abs(xin))))
dxin = step*eye(length(xin));
for k = 1:length(xin)
    df_fd(:,k) = (ipzt(obj,act,xin+dxin(:,k)) - ...
        ipzt(obj,act,xin-dxin(:,k)))/(2*step);
end


disp('Inverse pzt model error:');
disp(max(abs(df_fd(:)-df(:))));
% max(abs(df(:)))


end

function [V, dV] = ipzt(obj, act, xin)

nq = obj.getNumPositions();
nu = obj.getNumInputs();

t = xin(1);
x = xin(2:2*nq+1);
u = xin(2*nq+2:2*nq+1+nu);

act_joints = obj.getActuatedJoints();
qact = x(act_joints);     % convert back to SI

% % Drive Parameters
% Vb = 225;
% Vg = 0;
%
% % Geometry
% tact = 325e-3;          %m (actuator thickness)
% tpzt = 135e-3;          %m (pzt thickness)
% wn = 3.5;            %m (nominal width)
% wr = 1.5;               %m (width ratio)
% lext = 1;            %m (extension length)
% lact = 9.2;          %m (actuator length)
%
% % PZT internal Stress as a function of field (see N.T. Jafferis et al. 2015)
% f31_eff_b = 24e-3;         %Pa*m/V (in blocked case)
% f31_eff_f = 32e-3;         %Pa*m/V (in free case)
%
% % PZT Young's modulus as a function of strain (see N.T. Jafferis et al. 2015)
% Emin = 38.5e3;          %sigmoid min (Pa)
% Emax = 81e3;            %sigmoid max (Pa)
% Eave = 40e3;            % average modulus (Pa)
%
% % Material properties (see N.T. Jafferis et al. 2015)
% Ecf = 340e3;             % Pa (CF young's modulus)
%
% %  ----------- fillin in derived actuator properties -------------
% tcf = tact - 2*tpzt;            % cf thickness
% lr = lext/lact;                 % length ratio
% eps = (Vb -Vg)/tpzt;        % total field
%
% % geometric factor
% GF_num = 8*(1-wr)^3*(1 + 2*lr);
% GF_den = -6*(wr - 1)*(-3 + 4*lr*(wr - 1)+2*wr) + 3*(-2 + 2*lr*(wr-1)+wr)^2*log((2-wr)/wr);
% GF = GF_num/GF_den;
%
% % free deflection
% df_den = (1/3)*Eave*tpzt*(1.5*tcf^2 + 3*tcf*tpzt + 2*tpzt^2) + Ecf*tcf^3/12;
% df_num = f31_eff_f*tpzt*lact^2*eps*(tpzt+tcf);
% df = 0.25*(1 + 2*lr)*(df_num/df_den);

V = zeros(numel(u), 1);
dVdx = zeros(numel(u), numel(x));
dVdu = zeros(numel(u));

for i = 1:numel(qact)
    obj = act.dummy_bender(i); 
    Ct = (0.75/obj.ap.lact)*(obj.ap.f31_eff_b + (-obj.orien*qact(i)/obj.df)*...
        (obj.ap.f31_eff_f - obj.ap.f31_eff_b))*obj.ap.wn*(obj.ap.tpzt+obj.tcf)*obj.GF;
    Cb = (0.75/obj.ap.lact)*(obj.ap.f31_eff_b + (obj.orien*qact(i)/obj.df)*...
        (obj.ap.f31_eff_f - obj.ap.f31_eff_b))*obj.ap.wn*(obj.ap.tpzt+obj.tcf)*obj.GF;
    
    V(i) = (Ct*obj.dp.Vb + Cb*obj.dp.Vg - obj.orien*u(i))/(Ct+Cb);
    
    
    dCt_dq = (0.75/obj.ap.lact)*((-obj.orien/obj.df)*(obj.ap.f31_eff_f - obj.ap.f31_eff_b))*obj.ap.wn*(obj.ap.tpzt+obj.tcf)*obj.GF;
    dCb_dq = (0.75/obj.ap.lact)*((obj.orien/obj.df)*(obj.ap.f31_eff_f - obj.ap.f31_eff_b))*obj.ap.wn*(obj.ap.tpzt+obj.tcf)*obj.GF;
    
    d1dq = obj.dp.Vb*dCt_dq/(Ct + Cb);
    d2dq = obj.dp.Vg*dCb_dq/(Ct + Cb);
    d3dq = (u(i)*obj.orien*(dCt_dq + dCb_dq))/(Ct + Cb)^2;

    dVdx(i, act_joints(i)) = d1dq + d2dq + d3dq;    
    
    %d/du
    dVdu(i,i) = -obj.orien/(Ct + Cb);
end

dV = [zeros(numel(V), 1), dVdx, dVdu];
end
