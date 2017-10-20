function y = foamy_sensors(x,xdot)

%The Pixhawk Gazebo simulator uses a North-West-Up coordinate frame, so
%we've replicated that here

%Constants
RE = 6371000.0; %6353000.0 %Earth radius in meters
BASELAT = 42.467332*(pi/180); %47.397742; %Base station latitude
BASELON = -71.414699*(pi/180); %8.545594; %Base station longitude
BASEALT = 413.0; %488.0; %Base station altitude in meters

RHO = 1.225; %1.2754; %Density of air in kg/m^3 at sea level and 15C
%PRESS0 = 1013.2; %Sea level pressure in millibar

B0 = [0.194701; 0.050497; -0.480388]; %From NOAA for Concord flying field
%B0 = [0.21523; -0.00771; -0.42741]; %Taken from Pixhawk Gazebo code for Zurich in 10^5 nT, NWU

g = [0; 0; -9.81]; %NWU gravity field vector

x_rad = x(1)/RE;
y_rad = -x(2)/RE;
c = sqrt(x_rad*x_rad + y_rad+y_rad);
if c > 0
    sin_c = sin(c);
    cos_c = cos(c);
    lat_rad = asin(cos_c * sin(BASELAT) + (x_rad * sin_c * cos(BASELAT)) / c);
    lon_rad = (BASELON + atan2(y_rad*sin_c, c*cos(BASELAT)*cos_c - x_rad*sin(BASELAT)*sin_c));
else
    lat_rad = BASELAT;
    lon_rad = BASELON;
end

alt = BASEALT + x(3); %altitude in meters

%Velocity in m/s
vnorth = x(8);
veast = -x(9);
vdown = -x(10);

gps = [lat_rad; lon_rad; alt; vnorth; veast; vdown]; %GPS position and velocity
accel = qrotate(qconj(x(4:7)), xdot(8:10)-g); %Body frame acceleration
gyro = x(11:13); %body-frame angular velocity
mag = qrotate(qconj(x(4:7)), B0); %body-frame magnetic field vector

press = 0.0; %PRESS0 - .12*alt; %Atmospheric pressure in millibar
palt = x(3); %Pressure altitude in meters

v_body = qrotate(qconj(x(4:7)), x(8:10)); %Body frame velocity
pdyn = 0.01*0.5*RHO*(v_body(1)+v_body(3))*(v_body(1)+v_body(3)); %Dynamic pressure in millibar

y = [gps; accel; gyro; mag; press; pdyn; palt];

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
