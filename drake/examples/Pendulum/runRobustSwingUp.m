close all

p = PendulumPlant();
v = PendulumVisualizer();

N = 51;

[utraj1,xtraj1,z1,prog1] = p.swingUpDirtran(N);

%Playback
olsys1 = cascade(utraj1,p);
xol1 = olsys1.simulateODE(utraj1.tspan,[0 0]');
v.playback(xol1);

D = (1/.2^2); %This corresponds to +/-.1 uncertainty in mass (10%)
[utraj2,xtraj2,z2,prog2] = p.robustSwingUpTrajectory(N,D);

%Playback
olsys2 = cascade(utraj2,p);
xol2 = olsys2.simulateODE(utraj2.tspan,[0 0]');
v.playback(xol2);


%Closed-loop simulation
Q = [100 0; 0 10];
R = 1;
Qf = 1000*eye(2);

%p = p.setMass(1);
c1 = tvlqr(p,xtraj1,utraj1,Q,R,Qf);
p = p.setMass(1.2);
clsys1 = feedback(p,c1);
xcl1 = clsys1.simulate([utraj1.tspan(1) utraj1.tspan(2)+.5], [0 0]');
v.playback(xcl1);

p = p.setMass(1);
c2 = tvlqr(p,xtraj2,utraj2,Q,R,Qf);
p = p.setMass(1.2);
clsys2 = feedback(p,c2);
xcl2 = clsys2.simulate([utraj2.tspan(1) utraj2.tspan(2)+.5], [0 0]');
v.playback(xcl2);


% %Write movie files
% %v.playbackAVI(xcl1, 'swing1.avi');
% %v.playbackAVI(xcl2, 'swing2.avi');
% % setenv('PATH', [getenv('PATH') ':/usr/local/bin']);
% % setenv('PATH', [getenv('PATH') ':/Library/TeX/texbin']);
% % v.playbackSWF(xcl1, 'swing1.swf');
% % v.playbackSWF(xcl2, 'swing2.swf');
% 
% %Plots
% h1 = z1(prog1.h_inds);
% t1 = [0; cumsum(h1)];
% x1 = z1(prog1.x_inds);
% u1 = z1(prog1.u_inds);
% 
% h2 = z2(prog2.h_inds);
% t2 = [0; cumsum(h2)];
% x2 = z2(prog2.x_inds);
% u2 = z2(prog2.u_inds);
% 
% figure(2);
% subplot(2,1,1);
% plot(t1(1:end-1),u1);
% ylabel('u_{dirtran}');
% xlim([0 t2(end)]);
% ylim([-3.5 3.5]);
% subplot(2,1,2);
% plot(t2(1:end-1),u2);
% xlim([0 t2(end)]);
% ylim([-3.5 3.5]);
% ylabel('u_{robust}');
% 
% figure(3);
% subplot(2,1,1);
% plot(t1,x1(1,:));
% hold on
% plot(t1,x1(2,:));
% ylabel('x_{dirtran}');
% xlim([0 t2(end)]);
% l = legend('$\theta$', '$\dot{\theta}$');
% set(l,'Interpreter','latex')
% subplot(2,1,2);
% plot(t2,x2(1,:));
% hold on
% plot(t2,x2(2,:));
% ylabel('x_{robust}');
% xlim([0 t2(end)]);