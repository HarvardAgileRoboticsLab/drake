function [ output_args ] = foamy_pend_plotter(xtraj,utraj)


figure(10);
tspan1 = [0 xtraj.traj{1}.tspan(2)];
tspan2 = [xtraj.traj{2}.tspan(1) xtraj.traj{2}.tspan(2)];

ts1 = linspace(tspan1(1),tspan1(2));
ts2 = linspace(tspan2(1),tspan2(2));

xtraj1 = xtraj.traj{1}.trajs{2}.eval(ts1);
xtraj2 = xtraj.traj{2}.trajs{2}.eval(ts2);

%%Airplane

subplot(4,3,1);
plot(ts1,xtraj1(1,:));
hold on
plot(ts2,xtraj2(1,:));
title('Plane X');
hold off

subplot(4,3,2);
plot(ts1,xtraj1(2,:));
hold on
plot(ts2,xtraj2(2,:));
title('Plane Y');
hold off

subplot(4,3,3);
plot(ts1,xtraj1(3,:));
hold on
plot(ts2,xtraj2(3,:));
title('Plane Z');
hold off

subplot(4,3,4);
plot(ts1,xtraj1(15,:));
hold on
plot(ts2,xtraj2(15,:));
title('Plane xdot');
hold off

subplot(4,3,5);
plot(ts1,xtraj1(16,:));
hold on
plot(ts2,xtraj2(16,:));
title('Plane ydot');
hold off

subplot(4,3,6);
plot(ts1,xtraj1(17,:));
hold on
plot(ts2,xtraj2(17,:));
title('Plane zdot');
hold off

%%Hook
subplot(4,3,7);
plot(ts1,xtraj1(1,:));
hold on
plot(ts2,xtraj2(1,:));
title('Pend X');
hold off

subplot(4,3,8);
plot(ts1,xtraj1(2,:));
hold on
plot(ts2,xtraj2(2,:));
title('Pend Y');
hold off

subplot(4,3,9);
plot(ts1,xtraj1(3,:));
hold on
plot(ts2,xtraj2(3,:));
title('Plane Z');
hold off

subplot(4,3,10);
plot(ts1,xtraj1(15,:));
hold on
plot(ts2,xtraj2(15,:));
title('Pend xdot');
hold off

subplot(4,3,11);
plot(ts1,xtraj1(16,:));
hold on
plot(ts2,xtraj2(16,:));
title('Pend ydot');
hold off

subplot(4,3,12);
plot(ts1,xtraj1(17,:));
hold on
plot(ts2,xtraj2(17,:));
title('Pend zdot');
hold off

end


% 
%     r = x(1:3); %Lab-frame position vector of fuselage
%     q = x(4:7); %Quaternion rotation from body to lab frame for fuselage
%     R = qtoR(q);
%     r_p = x(8:10); %Lab-frame position vector of pendulum
%     q_p = x(11:14); %Quaternion rotation from body to lab frame for pendulum
%     R_p = qtoR(q_p);
%     v = x(15:17); %Lab-frame velocity vector of fuselage
%     w = x(18:20); %Body-frame angular velocity of fuselage
%     v_p = x(21:23); %Lab-frame velocity vector of pendulum
%     w_p = x(24:26); %Body-frame angular velocity of pendulum