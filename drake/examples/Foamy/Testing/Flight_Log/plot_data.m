function plot_data(veh,ts)
%PLOT_DATA Summary of this function goes here
%   Detailed explanation goes here
figure(10);

length = (ts(2)-ts(1))+1;
rpy_meas = zeros(length,3);
rpy_des = zeros(length,3);
for i = 1:length
    rpy_meas(i,:) = quat2rpy([veh.x_meas(i-1+ts(1),4),veh.x_meas(i-1+ts(1),5),veh.x_meas(i-1+ts(1),6),veh.x_meas(i-1+ts(1),7)]);
    rpy_des(i,:) = quat2rpy([veh.xd(i-1+ts(1),4),veh.xd(i-1+ts(1),5),veh.xd(i-1+ts(1),6),veh.xd(i-1+ts(1),7)]);
end


subplot(4,4,1);
plot(veh.times(ts(1):ts(2)),(veh.x_meas(ts(1):ts(2),1)+veh.x_init(ts(1):ts(2),1)));
hold on
plot(veh.times(ts(1):ts(2)),veh.xd(ts(1):ts(2),1));
title('X');
hold off

subplot(4,4,2);
plot(veh.times(ts(1):ts(2)),(veh.x_meas(ts(1):ts(2),2)+veh.x_init(ts(1):ts(2),2)));
hold on
plot(veh.times(ts(1):ts(2)),veh.xd(ts(1):ts(2),2));
title('Y');
hold off

subplot(4,4,3);
plot(veh.times(ts(1):ts(2)),(veh.x_meas(ts(1):ts(2),3)+veh.x_init(ts(1):ts(2),3)));
hold on
plot(veh.times(ts(1):ts(2)),veh.xd(ts(1):ts(2),3));
title('Z');
hold off

subplot(4,4,4);
 plot(veh.times(ts(1):ts(2)),rpy_meas(:,1));
 hold on
 plot(veh.times(ts(1):ts(2)),rpy_des(:,1)+2*pi);
%plot(veh.times(ts(1):ts(2)),veh.x_meas(ts(1):ts(2),6));
%hold on
%plot(veh.times(ts(1):ts(2)),veh.xd(ts(1):ts(2),6));
title('Roll');
hold off

subplot(4,4,5);
plot(veh.times(ts(1):ts(2)),rpy_meas(:,2));
hold on
plot(veh.times(ts(1):ts(2)),rpy_des(:,2));
title('Pitch');
hold off

subplot(4,4,6);
plot(veh.times(ts(1):ts(2)),rpy_meas(:,3));
hold on
plot(veh.times(ts(1):ts(2)),rpy_des(:,3));
title('Yaw');
hold off

%%Hook
subplot(4,4,7);
plot(veh.times(ts(1):ts(2)),veh.x_meas(ts(1):ts(2),8));
hold on
plot(veh.times(ts(1):ts(2)),veh.xd(ts(1):ts(2),8));
title('xdot');
hold off

subplot(4,4,8);
plot(veh.times(ts(1):ts(2)),veh.x_meas(ts(1):ts(2),9));
hold on
plot(veh.times(ts(1):ts(2)),veh.xd(ts(1):ts(2),9));
title('ydot');
hold off

subplot(4,4,9);
plot(veh.times(ts(1):ts(2)),veh.x_meas(ts(1):ts(2),10));
hold on
plot(veh.times(ts(1):ts(2)),veh.xd(ts(1):ts(2),10));
title('zdot');
hold off

subplot(4,4,10);
plot(veh.times(ts(1):ts(2)),veh.x_meas(ts(1):ts(2),11));
hold on
plot(veh.times(ts(1):ts(2)),veh.xd(ts(1):ts(2),11));
title('wx');
hold off

subplot(4,4,11);
plot(veh.times(ts(1):ts(2)),veh.x_meas(ts(1):ts(2),12));
hold on
plot(veh.times(ts(1):ts(2)),veh.xd(ts(1):ts(2),12));
title('wy');
hold off

subplot(4,4,12);
plot(veh.times(ts(1):ts(2)),veh.x_meas(ts(1):ts(2),13));
hold on
plot(veh.times(ts(1):ts(2)),veh.xd(ts(1):ts(2),13));
title('wz');
hold off

subplot(4,4,13);
plot(veh.times(ts(1):ts(2)),veh.u_comm(ts(1):ts(2),1));
title('U[0] Thr');


subplot(4,4,14);
plot(veh.times(ts(1):ts(2)),veh.u_comm(ts(1):ts(2),2));
title('U[1] Roll');


subplot(4,4,15);
plot(veh.times(ts(1):ts(2)),veh.u_comm(ts(1):ts(2),3));
title('U[2] Pitch');


subplot(4,4,16);
plot(veh.times(ts(1):ts(2)),veh.u_comm(ts(1):ts(2),4));
title('U[3] Yaw');



end

