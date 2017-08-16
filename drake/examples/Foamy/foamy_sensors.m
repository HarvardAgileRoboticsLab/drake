function y = foamy_sensors(x,u)

B0 = [-.0505; .1947; -.4807]; %ENU magnetic field vector in Gauss
g = [0; 0; -9.81]; %ENU gravity field vector

xdot = foamy_dynamics(0,x,u);

gps = [x(1:3); x(8:10)]; %Inertial position and velocity (East-North-Up)
accel = qrotate(qconj(x(4:7)), xdot(8:10)+g); %Body frame acceleration
gyro = x(11:13); %body-frame angular velocity
mag = qrotate(qconj(x(4:7)),B0); %body-frame magnetic field vector

y = [gps; accel; gyro; mag];

end

% ---------- Helper Functions ---------- %
function rrot = qrotate(q,r)
    %Rotate vector r by quaternion q
    rrot = r + 2*cross(q(2:4),(cross(q(2:4),r) + q(1)*r));
end

function qc = qconj(q)
    %Quaternion conjugate
    qc = [q(1); -q(2:4)];
end
