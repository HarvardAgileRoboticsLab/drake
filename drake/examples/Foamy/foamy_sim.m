function foamy_sim()

p = FoamyPlant();
[xtraj,utraj] = p.runDircol();

ts = .02;
Q = (1/.05)^2*eye(12);
Qn = (1/.01)^2*eye(12);
R = (1/10)^2*eye(4);

[x0,u0,K,t] = foamy_tvlqr(p,Q,R,Qn,xtraj,utraj,ts);

xtrue = zeros(13,length(xtraj)+1);
xtrue(:,1) = x0(:,1) + .05*randn(13,1);
xtrue(4:7,1) = xtrue(4:7,1)/norm(xtrue(4:7,1));
for k = 1:(length(t)-1);
    xtrue(:,k+1) = rkstep(xtrue(:,k),foamy_controller(xtrue(:,k),x0(:,k),u0(:,k),K(:,:,k)),ts);
end

xsim = PPTrajectory(foh(t,xtrue));
xsim = xsim.setOutputFrame(p.getStateFrame());
v = FoamyVisualizer(p);
v.playback(xsim);

end

function xn = rkstep(x,u,dt)
    %Single step of the 3rd order Runge-Kutta method
    xdot1 = foamy_dynamics_mex(0,x,u);
    xdot2 = foamy_dynamics_mex(0,x+.5*dt*xdot1,u);
    xdot3 = foamy_dynamics_mex(0,x-dt*xdot1+2*dt*xdot2,u);
    xn = x + (dt/6)*(xdot1 + 4*xdot2 + xdot3);
end