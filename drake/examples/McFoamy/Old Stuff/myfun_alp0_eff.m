function F = myfun_alp0_eff(x,Kp,Kv,delCL)
alp = 0;
alp0_eff = x;

a_eff = alp - alp0_eff;

F = Kp*sin(a_eff)*(cos(a_eff))^2 + Kv*abs(sin(a_eff))*sin(a_eff)*cos(a_eff) - delCL;