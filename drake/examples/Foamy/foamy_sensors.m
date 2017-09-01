function y = foamy_sensors(x,xdot)

%Constants
RE = 6371000.0; %Earth radius in meters
BASELAT = 47.3977420; %42.467332; %Base station latitude in degrees
COSLAT = 0.676904978400579; %0.737662414272813; %cosine of base station latitude
BASELON = 8.5455937; %-71.414699; %Base station longitude in degrees
BASEALT = 488.11; %413.0; %Base station altitude in meters

RHO = 1.225; %Density of air in kg/m^3 at sea level and 15C
PRESS0 = 1013.2; %Sea level pressure in millibar

B0 = [.008665; .21489; -.428427]; %[-.0505; .1947; -.4807]; %ENU magnetic field vector in Gauss
g = [0; 0; -9.81]; %ENU gravity field vector

lat = BASELAT + (180.0/pi)*(x(2)/RE); %lat in degrees (linearized)
lon = BASELON + (180.0/pi)*(x(1)/(RE*COSLAT)); %lon in degrees (linearized)
alt = BASEALT + x(3); %altitude in meters

%Velocity in m/s
vnorth = x(9);
veast = x(8);
vdown = -x(10);

gps = [lat; lon; alt; vnorth; veast; vdown]; %GPS position and velocity
accel = qrotate(qconj(x(4:7)), xdot(8:10)-g); %Body frame acceleration
gyro = x(11:13); %body-frame angular velocity
mag = qrotate(qconj(x(4:7)),B0); %body-frame magnetic field vector

press = PRESS0 - .12*alt; %Atmospheric pressure in millibar
pdyn = 0.01*0.5*RHO*(x(8:10)'*x(8:10)); %Dynamic pressure in millibar

y = [gps; accel; gyro; mag; press; pdyn];

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
