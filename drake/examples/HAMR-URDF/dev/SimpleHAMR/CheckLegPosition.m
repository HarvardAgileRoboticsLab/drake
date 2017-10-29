clear; close all; clc;

%%
sl_urdf = fullfile(getDrakePath, 'examples', 'HAMR-URDF', 'dev', 'SimpleHAMR', 'urdf', 'FL_scaled.urdf');

% options
% options.terrain = RigidBodyFlatTerrain(); 
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = false;
% options.collision = true; 

SL = RigidBodyManipulator(sl_urdf, options);
SL = SL.setGravity([0; 0; -9.81e-3]);
SL = compile(SL);

nq = SL.getNumPositions();
nv = SL.getNumVelocities();
nx = nq + nv;
nu = SL.getNumInputs();
nl = SL.getNumStateConstraints();

pfSL = [0 7.58, -11.35]'; % position of foot in local frame    

v = SL.constructVisualizer();
x0 = zeros(nx,1); %x0(3) = 13; 
v.inspector(x0);

fkoptSL.base_or_frame_id = SL.findLinkId('Chassis');
kinsoliSL = SL.doKinematics(x0(1:nq), x0(1:nq)*0);
xfootiSL = SL.forwardKin(kinsoliSL, SL.findLinkId('FLL4'), pfSL , fkoptSL); 

%%
hamr_urdf = fullfile(getDrakePath, 'examples', 'HAMR-URDF', 'dev', 'SimpleHAMR', ...
    'urdf', 'HAMRSimple_scaled.urdf');

% options
% options.terrain = RigidBodyFlatTerrain(); 
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = true;


hamr = HAMRSimpleRBM(hamr_urdf,options);
x0 = hamr.getInitialState(); 
% v = hamr.constructVisualizer();
% v.inspector(x0); 

nq = hamr.getNumPositions();
nv = hamr.getNumVelocities();
nc = hamr.getNumContactPairs();

pf = [0 0 -14.988382167532292]'; % position of foot in local frame    

fkopt.base_or_frame_id = hamr.findLinkId('Chassis');

kinsoli = hamr.doKinematics(x0(1:nq), 0*x0(1:nq));
xfootiHAMR = hamr.forwardKin(kinsoli, hamr.findLinkId('FL2'), pf, fkopt); 

xfootiSL - xfootiHAMR