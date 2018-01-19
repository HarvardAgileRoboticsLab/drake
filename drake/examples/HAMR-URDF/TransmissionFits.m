clear; clc; close all; 

npts = 11;
xlim = [-3, 3];
zlim = [-1.5, 3];

% yleg: 3rd dim (x, y, z) and 4th dim (FL, RL, FR, RR);
% H: 3rd dim (FL, RL, FR, RR);
[xx, zz, qact, H] = gen_trans_data(npts, xlim, zlim);

save('ActuatorToLegVelocity', 'xx', 'zz', 'H')
