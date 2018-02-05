clear; clc; %close all;

% actuator
% actuator = H6_V3_actuator(325e-6); 
actuator = H6_V3_actuator_scaled(335e-3); 

% extract actuator geometry 
lext = actuator.lext;          % extension length (mm)
L = actuator.l;                % pzt length (mm)
wn = actuator.wn;              % nominal width (mm)
wr = actuator.wr;
lr = lext/L; 

% PZT width
tpzt = actuator.tpzt;          %pzt-thickness (mm)  

f31_eff_b = actuator.f31_eff_b;         %Pa*mm/V (in blocked case)
f31_eff_f = actuator.f31_eff_f;         %Pa*mm/V (in free case)

% PZT internal Stress as a function of field 
% f31_min = actuator.f31_min;         %Pa*m/V (sigmoid min)
% f31_max = actuator.f31_max;         %Pa*m/V (sigmoid max)
% field0 = actuator.field0;           %V/m (switching field)
% c = actuator.c;                     %sigmoid width
% d = actuator.d;                     %m/V (linear dependence of stress on field at high fields)

% PZT Young's modulus as a function of strain
Eave_eff = actuator.Eave_eff;
% Emin = actuator.Emin;          %sigmoid min (Pa)
% Emax = actuator.Emax;          %sigmoid max (Pa)
% eps0 = actuator.eps0;          %switching strain
% a = actuator.a;                %sigmoid width

% linear depedence of stress on strain 
% b = actuator.b; 

% CF propreties
tcf = actuator.tcf;       % carbon fiber thickness (m);
Ecf = actuator.Ecf;       % young's modulus (Pa)

% geometric factor
GF_num = 8*(1-wr)^3*(1 + 2*lr);
GF_den = -6*(wr - 1)*(-3 + 4*lr*(wr - 1)+2*wr) + ...
    3*(-2 + 2*lr*(wr-1)+wr)^2*log((2-wr)/wr);
GF = GF_num/GF_den;

Vb = 225; 
field = Vb/tpzt; 



Fb = 0.75*f31_eff_b*wn*tpzt*(tpzt+tcf)*GF*field/L
df_num = f31_eff_f*tpzt*(tpzt+tcf)*L^2*field; 
df_den = (1/3)*Eave_eff*tpzt*(1.5*tcf^2 + 3*tcf*tpzt + 2*tpzt^2) + Ecf*tcf^3/12; 
df_pp = 0.5*(df_num/df_den)*(1 + 2*lr)

k = 2*Fb/df_pp
k2 = 3*(f31_eff_b/f31_eff_f)*(wn/L^3)*GF*(df_den/(1+2*lr))
% 
% q = linspace(-1, 1, 1000); 
% % tcf = linspace(0, 150e-3, 100); 
% % [q, tcf] = meshgrid(q, tcf); 
% eps = q*((tpzt+tcf)/(2*(1 + 2*lr)*L^2)); 
% sig1 = eps*Emin - ((Emax-Emin)/a)*(log(1 + exp(a*(eps0-eps)))-log(1+exp(a*eps0))); 
% sig2 = -eps*Emin - ((Emax-Emin)/a)*(log(1 + exp(a*(eps0+eps)))-log(1+exp(a*eps0))); 
% sig = sig1 - sig2; 
% % sig_ave_2 = eps*Emin-((Emax-Emin)/(2*a))*log((1 + exp(a*(eps0-eps)) + exp(a*(eps0+eps)) + exp(2*a*eps0))/...
% %     (1 + 2*exp(a*eps0) + exp(2*a*eps0))); 
% 
% % eps2 = -abs(q)*((tpzt+tcf)/(2*(1 + 2*lr)*L^2)); 
% eps_tot = -abs(q)*((tpzt+tcf)/(2*(1 + 2*lr)*L^2)); 
% 
% % delta = 1e-8; 
% % eps_soft = -(sqrt(q.^2 + delta^2))*(tpzt+tcf)/(2*(1 + 2*lr)*L^2); 
% figure(1); hold on; 
% plot(q, sig1);
% plot(q, sig2);
% plot(q, sig);
% % plot(q, sig_ave_2, '--');
% % % plot(q, eps_soft, '--'); 
% % legend('True', 'Fit'); 
% % 
% % 
% % % 
% % % 
% % bs_term = 1./(2*a*eps_tot).* (log(1 + exp(a*(eps0 - eps_tot))) - log(1 + exp(a*eps0))); 
% % 
% % 
% % 
% % a0 = -0.1659;
% a1 = 0.1561;
% w = 2.233; 
% 
% bs_fit = a0 + a1*cos(q*w);
% % 
% [pft, ~, mu] = polyfit(q, bs_term, 8);
% bs_fit = polyval(pft, (q-mu(1))/mu(2)); 
% % 
% % surf(q, tcf, bs_term); 
% % xlabel('Deflection')
% % ylabel('Thickness')

% c = exp(a*(eps0-eps));
% E = Emin + (Emax - Emin)*(c./(1+c)); 

% c2 = exp(a*(eps0-eps2));
% E2  = Emin + (Emax - Emin)*(c2./(1+c2)); 
% Eave =  Emin - ((Emax - Emin)./(a*eps_tot)).*...
%     (log(1 + exp(a*(eps0 - eps_tot))) - log(1 + exp(a*eps0))); 
% 
% % electric field
% % Eave_fit = Emin - ((Emax - Emin))*bs_fit; 
% K2 = 0.8*3*(wn/L^3/(1+2*lr))*GF*(Ecf*tcf.^3/12).*q + ...
%     0.8*3*(wn/L)*GF*(((1/3)*sig*tpzt*(1.5*tcf^2 + 3*tcf*tpzt + 2*tpzt^2))/(tpzt+tcf)); 
% 
% 
%     %(0.85*3*(wn/L/(tpzt+tcf))*GF*((1/3)*sig_ave*tpzt*(1.5*tcf^2 + 3*tcf*tpzt + 2*tpzt^2)))./q; 
% 
% K = 0.85*3*(wn/L^3)*GF*(((1/3)*Eave*tpzt*(1.5*tcf^2 + 3*tcf*tpzt + 2*tpzt^2)+ ...
%     Ecf*tcf.^3/12)/(1 + 2*lr));             % actuator stiffness
% % K2 = 0.85*3*(wn/L^3)*GF*(((1/3)*E2*tpzt*(1.5*tcf^2 + 3*tcf*tpzt + 2*tpzt^2)+ ...
% %     Ecf*tcf.^3/12)/(1 + 2*lr));             % actuator stiffness
% % Kfit = 0.85*3*(wn/L^3)*GF*(((1/3)*Eave_fit*tpzt*(1.5*tcf^2 + 3*tcf*tpzt + 2*tpzt^2)+ ...
% %     Ecf*tcf^3/12)/(1 + 2*lr));     
% % % 
% % 
% figure(3); clf; hold on; 
% % subplot(1,3,1); 
% % plot(q, bs_term); hold on; 
% % plot(q, bs_fit); hold on; 
% % title('Bs Term')
% % plot(eps_soft, bs_fit, '--'); 
% subplot(1,3,2); 
% % plot(q, E); hold on; 
% % plot(-q, E); hold on; 
% % plot(q, E2, '--'); 
% plot(q, Eave.*eps); hold on; 
% plot(q, sig)
% % plot(q, Eave_fit); 
% xlabel('Tip Deflection (mm)')
% % plot(eps_soft, Eave_fit, '--'); 
% title('Stress')
% subplot(1,3,3)
% plot(q, K); hold on; 
% plot(q, K2./q, '--'); hold on; 
% % plot(q, Kfit); hold on; 
% title('Force')
% xlabel('Tip Deflection (mm)')
% % plot(eps_soft, Kfit); 
% 
% % figure(2); clf; hold on; 
% % plot(q, eps_tot); 