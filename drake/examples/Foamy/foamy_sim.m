function foamy_sim()

p = FoamyPlant();
[xtraj,utraj] = p.runDircol();

ts = .02;
Q = (1/.05)^2*eye(12);
Qn = (1/.01)^2*eye(12);
R = (1/.25)^2*eye(4);

[x0,u0,K,t] = foamy_tvlqr(p,Q,R,Qn,xtraj,utraj,ts);

Qkf = (.001)^2*eye(12);
Rkf = (.05)^2*eye(17);
Phat = (.1)^2*eye(12);

xtrue = zeros(13,length(t));
xtrue(:,1) = x0(:,1);
xtrue(3,1) = xtrue(3,1)-1;
xtrue(2,1) = xtrue(2,1)+3;
xtrue(4:7,1) = xtrue(4:7,1)/norm(xtrue(4:7,1));
xhat = x0(:,1);
for k = 1:(length(t)-1);
    u = foamy_controller(xhat,x0(:,k),u0(:,k),K(:,:,k));
    [xtrue(:,k+1), xdot] = foamy_rk4(xtrue(:,k),u,ts);
    y = foamy_sensors(xtrue(:,k+1),xdot);
    [xhat, Phat] = foamy_ukf(xhat,y,u,Phat,Qkf,Rkf,ts);
end

xsim = PPTrajectory(foh(t,xtrue));
xsim = xsim.setOutputFrame(p.getStateFrame());
v = FoamyVisualizer(p);
v.playback(xsim,struct('slider',true));

end
