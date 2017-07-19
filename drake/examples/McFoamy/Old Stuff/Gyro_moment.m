function [Gyro_moment] = Gyro_moment(Ang_vel,B_rate,J)
%GYRO_MOMENT Summary of this function goes here

GM1 = 0;
GM2 = B_rate(3)*Ang_vel*J*(-1);
GM3 = B_rate(2)*Ang_vel*J;

Gyro_moment = [GM1;GM2;GM3];

end

