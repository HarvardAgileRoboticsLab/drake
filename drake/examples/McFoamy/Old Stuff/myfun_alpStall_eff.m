function F = myfun_alpStall_eff(x,Kp,Kv,alp0_eff,CLmax)
alp = x;

a_eff = alp - alp0_eff;
F = Kp*sin(a_eff)*(cos(a_eff))^2 + Kv*abs(sin(a_eff))*sin(a_eff)*cos(a_eff) - CLmax;