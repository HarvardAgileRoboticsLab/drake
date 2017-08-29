function xdot = foamy_dynamics(t,x,u)

    %TODO: Add body aerodynamic forces

    %State vector:
    %r = x(1:3); %Lab-frame position vector
    q = x(4:7); %Quaternion rotation from body to lab frame
    v = x(8:10); %Lab-frame velocity vector
    w = x(11:13); %Body-frame angular velocity

    %Control input:
    thr = u(1); %Throttle command (0 to 1 for Pixhawk)
    ail = u(2); %Aileron command (-1 to 1 for Pixhawk)
    elev = u(3); %Elevator command (-1 to 1 for Pixhawk)
    rud = u(4); %Rudder command (-1 to 1 for Pixhawk)

    %Note that body coordinate frame is:
    % x: points forward out nose
    % y: points out right wing tip
    % z: points down

    % ---------- Input Checks ---------- %
%     q = q/sqrt(q'*q); %make sure quaternion is normalized
%     thr = min(1, max(0, thr));
%     ail = min(1, max(-1, ail));
%     elev = min(1, max(-1, elev));
%     rud = min(1, max(-1, rud));

    % ---------- Model Parameters ---------- %
    p = foamy_parameters; %load model parameters

    % ---------- Map Control Inputs to Angles ---------- %
    delta_ail = (ail-p.trim_ail)*p.g_ail;
    delta_elev = (elev-p.trim_elev)*p.g_elev;
    delta_rud = (rud-p.trim_rud)*p.g_rud;

    % ---------- Aerodynamic Forces (body frame) ---------- %

    v_body = qrotate(qconj(q),v); %body-frame velocity
    v_propwash = propwash(thr);
    v_rout = v_body + cross(w,[0; p.r_ail; 0]);
    v_lout = v_body + cross(w,[0; -p.r_ail; 0]);
    v_rin = v_body + cross(w,[0; p.l_in; 0]) + v_propwash;
    v_lin = v_body + cross(w,[0; -p.l_in; 0]) + v_propwash;
    v_elev = v_body + cross(w,[-p.r_elev; 0; 0]) + v_propwash;
    v_rud = v_body + cross(w,[-p.r_rud; 0; -p.z_rud]) + v_propwash;

    % --- Outboard Wing Sections --- %
    a_rout = alpha(v_rout);
    a_lout = alpha(v_lout);
    a_eff_rout = a_rout + p.ep_ail*delta_ail; %effective angle of attack
    a_eff_lout = a_lout - p.ep_ail*delta_ail; %effective angle of attack

    F_rout = -p_dyn(v_rout)*.5*p.S_out*[Cd_wing(a_eff_rout); 0; Cl_wing(a_eff_rout)];
    F_lout = -p_dyn(v_lout)*.5*p.S_out*[Cd_wing(a_eff_lout); 0; Cl_wing(a_eff_lout)];

    F_rout = arotate(a_rout,F_rout); %rotate to body frame
    F_lout = arotate(a_lout,F_lout); %rotate to body frame

    % --- Inboard Wing Sections (Includes Propwash) --- %
    a_rin = alpha(v_rin);
    a_lin = alpha(v_lin);
    a_eff_rin = a_rin + p.ep_ail*delta_ail; %effective angle of attack
    a_eff_lin = a_lin - p.ep_ail*delta_ail; %effective angle of attack

    F_rin = -p_dyn(v_rin)*.5*p.S_in*[Cd_wing(a_eff_rin); 0; Cl_wing(a_eff_rin)];
    F_lin = -p_dyn(v_lin)*.5*p.S_in*[Cd_wing(a_eff_lin); 0; Cl_wing(a_eff_lin)];

    F_rin = arotate(a_rin,F_rin); %rotate to body frame
    F_lin = arotate(a_lin,F_lin); %rotate to body frame

    % --- Elevator --- %
    a_elev = alpha(v_elev);
    a_eff_elev = a_elev + p.ep_elev*delta_elev; %effective angle of attack

    F_elev = -p_dyn(v_elev)*p.S_elev*[Cd_elev(a_eff_elev); 0; Cl_elev(a_eff_elev)];

    F_elev = arotate(a_elev,F_elev); %rotate to body frame

    % --- Rudder --- %
    a_rud = beta(v_rud);
    a_eff_rud = a_rud + p.ep_rud*delta_rud; %effective angle of attack

    F_rud = -p_dyn(v_rud)*p.S_rud*[Cd_rud(a_eff_rud); Cl_rud(a_eff_rud); 0];

    F_rud = brotate(a_rud,F_rud); %rotate to body frame

    % --- Propeller --- %
    F_thr = [thr*p.g_thr; 0; 0];
    n_prop = sqrt(F_thr(1)/(p.rho*p.Ct*(p.D_prop^4))); %rotation speed in Hz
    w_prop = [2*pi*n_prop 0 0]';
    q_prop = p.rho*p.Cq*(n_prop^2)*(p.D_prop^5);
    T_prop = [-q_prop 0 0]';

    % ---------- Aerodynamic Torques (body frame) ---------- %

    T_rout = cross([0; p.r_ail; 0],F_rout);
    T_lout = cross([0; -p.r_ail; 0],F_lout);

    T_rin = cross([0; p.l_in; 0],F_rin);
    T_lin = cross([0; -p.l_in; 0],F_lin);

    T_elev = cross([-p.r_elev; 0; 0],F_elev);

    T_rud = cross([-p.r_rud; 0; -p.z_rud],F_rud);

    % ---------- Add Everything Together ---------- %

    F_aero = F_rout + F_lout + F_rin + F_lin + F_elev + F_rud + F_thr;
    F = qrotate(q,F_aero) - [0; 0; p.m*p.g];

    T = T_rout + T_lout + T_rin + T_lin + T_elev + T_rud + T_prop;

    xdot = [v;
            .5*[-q(2:4)'*w; q(1)*w + cross(q(2:4),w)];
            F/p.m;
            p.Jinv*(T - cross(w,(p.J*w + p.Jprop*w_prop)))];
end

function a = alpha(v)
    %Angle of attack
    a = atan2(v(3),v(1));
end

function b = beta(v)
    %Sideslip angle
    b = atan2(v(2),v(1));
end

function rrot = qrotate(q,r)
    %Rotate vector r by quaternion q
    rrot = r + 2*cross(q(2:4),(cross(q(2:4),r) + q(1)*r));
end

function qc = qconj(q)
    %Quaternion conjugate
    qc = [q(1); -q(2:4)];
end

function rrot = arotate(a,r)
    %Rotate by angle of attack
    rrot = [cos(a) 0  -sin(a);
              0    1    0;
            sin(a) 0  cos(a)]*r;
end

function rrot = brotate(b,r)
    %Rotate by sideslip angle
    rrot = [cos(b) -sin(b) 0;
            sin(b)  cos(b) 0;
              0       0    1]*r;
end

function v = propwash(thr)
    %Propwash wind speed (body frame)
    %From Bernoulli's equation
    
    p = foamy_parameters;
    
    v = [sqrt(2*thr*p.g_thr/(p.rho*p.A_prop)); 0; 0];
end

function pd = p_dyn(v)
    %Dynamic pressure
    
    p = foamy_parameters; %load model parameters
    
    pd = .5*p.rho*(v'*v);
end

function cl = Cl(a)

    a = min(pi/2, max(-pi/2, a));
    p = foamy_parameters;
    cl = polyval(p.Clcoef, a);
    
end

function cd = Cd(a)
    
    a = min(pi/2, max(-pi/2, a));
    p = foamy_parameters;
    cd = polyval(p.Cdcoef, a);

end

function cl = Cl_wing(a)
    %Lift coefficient (alpha in radians)
    
    cl = Cl(a);
%     a = min(pi/2, max(-pi/2, a));
%     
%     p = foamy_parameters; %load model parameters
%     
%     cl = (pi/2)*p.Ra*a; %flat plate theory for thin finite-length wing
end

function cd = Cd_wing(a)
    %Drag coefficient (alpha in radians)
    %Induced drag for a tapered finite wing
    %From Phillips P.55
    
    cd = Cd(a);
%     a = min(pi/2, max(-pi/2, a));
%     
%     p = foamy_parameters; %load model parameters
%     
%     cd = (pi/4)*p.Ra*a^2;
end

function cl = Cl_elev(a)
    %Lift coefficient (alpha in radians)
    
%     cl = Cl(a);
    
    a = min(pi/2, max(-pi/2, a));
    
    p = foamy_parameters; %load model parameters
    
    cl = (pi/2)*p.Ra_elev*a;
end

function cd = Cd_elev(a)
    %Drag coefficient (alpha in radians)
    %Induced drag for a tapered finite wing
    %From Phillips P.55

%     cd = Cd(a);
    
    a = min(pi/2, max(-pi/2, a));
    
    p = foamy_parameters; %load model parameters

    cd = (pi/4)*p.Ra_elev*a^2;
end

function cl = Cl_rud(a)
    %Lift coefficient (alpha in radians)
    
%     cl = Cl(a);
    
    a = min(pi/2, max(-pi/2, a));
    
    p = foamy_parameters; %load model parameters
    
    cl = (pi/2)*p.Ra_rud*a;
end

function cd = Cd_rud(a)
    %Drag coefficient (alpha in radians)
    %Induced drag for a tapered finite wing
    %From Phillips P.55

%     cd = Cd(a);
    
    a = min(pi/2, max(-pi/2, a));
    
    p = foamy_parameters; %load model parameters
    
    cd = (pi/4)*p.Ra_rud*a^2;
end