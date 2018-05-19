clear; clc; close all;
%
%% Build robot

name = 'HAMR_scaledV2_TYM';
urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', [name, '.urdf']);

options.terrain = RigidBodyFlatTerrain();
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = true;
options.collision = true;
options.dt = 1;

% Build robot + visualizer
hamr = HamrTSRBM(urdf, options);
nq = hamr.getNumPositions();
nv = hamr.getNumVelocities();
nu = hamr.getNumInputs();
v = hamr.constructVisualizer();
dim = 3; 
nc = 4; 

tsim = 200; 
x0 = hamr.getInitialState(); x0(3) = x0(3); 
q0 = x0(1:nq); 
utraj = PPTrajectory(zoh([0, tsim], zeros(nu,2))); 
utraj = utraj.setOutputFrame(hamr.getInputFrame()); 
disp('Simulating...')
xtraj = hamr.simulate([0, tsim], x0); 
disp('Done')

xx = xtraj.eval(xtraj.getBreaks());
xf = xx(:,end); 

%% Playback 
xtraj_scaled = DTTrajectory(xtraj.getBreaks()*1e-3, xtraj.eval(xtraj.getBreaks()));
xtraj_scaled = xtraj_scaled.setOutputFrame(xtraj.getOutputFrame());
options.slider = true;
v.playback(xtraj_scaled, options);

%% Marker Positions from Solidworks

BFR_MDL= [19.831, -9.790, 9.610];
BML_MDL = [0, 9.790, 4.721];
BRR_MDL = [-19.831, -9.790, 9.610];

FOOT_MDL = [0.06, 8.565, -6.322;
-0.06, 8.565, -6.322;
0.06, -8.565, -6.322;
-0.06, -8.565, -6.322];

legs = {'FLL4', 'RLL4', 'FRL4', 'RRL4'};

kinsol_AIR = hamr.doKinematics(q0, 0*q0);
body_marker_positions_MDL_AIR = hamr.forwardKin(kinsol_AIR, hamr.findLinkId('Chassis'), [BFR_MDL', BML_MDL', BRR_MDL']);
foot_marker_positions_MDL_AIR = zeros(dim, nc);

fkopt.base_or_frame_id = hamr.findLinkId('Chassis');
for i = 1:nc        
        foot_marker_positions_MDL_AIR(:,i)  = hamr.forwardKin(kinsol_AIR, hamr.findLinkId(legs{i}), FOOT_MDL(i,:)', fkopt);
end

kinsol_GND = hamr.doKinematics(xf(1:nq), xf(nq+(1:nv)));
body_marker_positions_MDL_GND = hamr.forwardKin(kinsol_GND, hamr.findLinkId('Chassis'), [BFR_MDL', BML_MDL', BRR_MDL']);
foot_marker_positions_MDL_GND = zeros(dim, nc);

for i = 1:nc        
        foot_marker_positions_MDL_GND(:,i)  = hamr.forwardKin(kinsol_GND, hamr.findLinkId(legs{i}), FOOT_MDL(i,:)', fkopt);
end


marker_positions_MDL_AIR = [body_marker_positions_MDL_AIR, foot_marker_positions_MDL_AIR]; 
marker_positions_MDL_AIR0 = bsxfun(@minus, marker_positions_MDL_AIR, mean(marker_positions_MDL_AIR,2));

marker_positions_MDL_GND = [body_marker_positions_MDL_GND, foot_marker_positions_MDL_GND]; 
marker_positions_MDL_GND0 = bsxfun(@minus, marker_positions_MDL_GND, mean(marker_positions_MDL_GND,2));

delta_model = marker_positions_MDL_GND - marker_positions_MDL_AIR

%% Marker Positions from Vicon (AIR)
% 5/16/2018
% FR_JBRN_AIR = [561.9, 358.8, 318.7]/10;
% ML_JBRN_AIR = [380.1, 538.6, 275.9]/10;
% RR_JBRN_AIR = [187.8, 345.0, 318.2]/10;
% FL_AIR = [500.9, 640.0, 138.4]/10;
% RL_AIR = [239.1, 629.1, 134.3]/10;
% FR_AIR = [523.1, 272.6, 124.3]/10;
% RR_AIR = [261.9, 266.1, 123.1]/10;

% 5/17/2018
FR_JBRN_AIR = [422.0, 79.9, 320.5]/10;
ML_JBRN_AIR = [251.2, 269.9, 276.6]/10;
RR_JBRN_AIR = [47.7, 87.6, 316.4]/10;
FL_AIR = [380.3, 364.2, 140.6]/10;
RL_AIR = [116.0, 368.5, 132.8]/10;
FR_AIR = [380.7, -3.4, 125.2]/10;
RR_AIR = [118.6, 4.3, 121.5]/10;

marker_positions_JBRN_AIR = [FR_JBRN_AIR', ML_JBRN_AIR', RR_JBRN_AIR', FL_AIR', RL_AIR', FR_AIR', RR_AIR'];
marker_positions_JBRN_AIR0 = bsxfun(@minus, marker_positions_JBRN_AIR, mean(marker_positions_JBRN_AIR,2));


%% Marker Positons from Vicon (GND)
% FR_JBRN_GND = [466.9, 404.5, 224.4]/10;
% ML_JBRN_GND = [283.1, 581.7, 179.3]/10;
% RR_JBRN_GND = [92.9, 388.1, 230.3]/10;
% FL_GND = [400.4, 683.5, 41.1]/10;
% RL_GND = [141.0, 672.6, 43.9]/10;
% FR_GND = [425.2, 307.0, 39.8]/10;
% RR_GND = [164.5, 300.5, 40.4]/10;

% 5/17/2018
FR_JBRN_GND = [363.2, 206.2, 220.3]/10;
ML_JBRN_GND = [193.0, 396.7, 175.9]/10;
RR_JBRN_GND = [-10.9, 218.0, 228.3]/10;
FL_GND = [316.1, 489.7, 37.7]/10;
RL_GND = [55.9, 497.4, 40.3]/10;
FR_GND = [313.7, 112.2, 35.6]/10;
RR_GND = [52.4, 123.6, 39.2]/10;

marker_positions_JBRN_GND = [FR_JBRN_GND', ML_JBRN_GND', RR_JBRN_GND', FL_GND', RL_GND', FR_GND', RR_GND'];
marker_positions_JBRN_GND0 = bsxfun(@minus, marker_positions_JBRN_GND, mean(marker_positions_JBRN_GND,2));

delta_jbrn = marker_positions_JBRN_GND0 - marker_positions_JBRN_AIR0


figure(1); clf; 
subplot(1,2,1); hold on; view(3);
plot3(marker_positions_MDL_AIR0(1,:),    marker_positions_MDL_AIR0(2,:), ...
    marker_positions_MDL_AIR0(3,:), 'k*')
plot3(marker_positions_MDL_GND0(1,:),    marker_positions_MDL_GND0(2,:), ...
    marker_positions_MDL_GND0(3,:), 'r*')
axis equal; 
subplot(1,2,2); hold on; view(3);
plot3(marker_positions_JBRN_AIR0(1,:),    marker_positions_JBRN_AIR0(2,:), ...
    marker_positions_JBRN_AIR0(3,:), 'k*')
plot3(marker_positions_JBRN_GND0(1,:),    marker_positions_JBRN_GND0(2,:), ...
    marker_positions_JBRN_GND0(3,:), 'r*')
axis equal;

