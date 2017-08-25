function y = foamy_sensors(x,u)

%Constants
RE = 6371000.0; %Earth radius in meters
BASELAT = 42.467332; %Base station latitude in degrees
COSLAT = 0.737662414272813; %cosine of base station latitude
BASELON = -71.414699; %Base station longitude in degrees
BASEALT = 413.0; %Base station altitude in meters

B0 = [-.0505; .1947; -.4807]; %ENU magnetic field vector in Gauss
g = [0; 0; -9.81]; %ENU gravity field vector

xdot = foamy_dynamics(0,x,u);

lat = BASELAT + (180.0/pi)*(x(2)/RE); %lat in degrees (linearized)
lon = BASELON + (180.0/pi)*(x(1)/(RE*COSLAT)); %lon in degrees (linearized)
alt = BASEALT + x(3); %altitude in meters

vlat = x(9);
vlon = x(8);
valt = x(10);

gps = [lat; lon; alt; vlat; vlon; valt]; %GPS position and velocity
accel = qrotate(qconj(x(4:7)), xdot(8:10)-g); %Body frame acceleration
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
