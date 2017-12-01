clear; close all; clc;

%%
hamr1_urdf = fullfile(getDrakePath, 'examples', 'HAMR-URDF', 'dev', 'SimpleHAMR', ...
    'urdf', 'HAMRSimple_scaled.urdf');
% options
% options.terrain = RigidBodyFlatTerrain(); 
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = false;
% options.collision = true; 

hamr1 = RigidBodyManipulator(hamr1_urdf, options);
hamr1 = hamr1.setGravity([0; 0; -9.81e-3]);
hamr1 = compile(hamr1);

nq = hamr1.getNumPositions();
nv = hamr1.getNumVelocities();
nx = nq + nv;
nu = hamr1.getNumInputs();
nl = hamr1.getNumStateConstraints();

pfSL = [0 7.58, -11.35]'; % position of foot in local frame    

v = hamr1.constructVisualizer();
x0 = zeros(nx,1); %x0(3) = 13; 
v.inspector(x0);

fkoptSL.base_or_frame_id = hamr1.findLinkId('Chassis');
kinsoliSL = hamr1.doKinematics(x0(1:nq), x0(1:nq)*0);
% xfootiSL = hamr1.forwardKin(kinsoliSL, hamr1.findLinkId('FLL4'), pfSL , fkoptSL); 

hamr1.getMass()

%%
hamr_urdf = fullfile(getDrakePath, 'examples', 'HAMR-URDF', 'urdf', 'HAMR_scaledV2.urdf');

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
% xfootiHAMR = hamr.forwardKin(kinsoli, hamr.findLinkId('FL2'), pf, fkopt); 

hamr.getMass
% xfootiSL - xfootiHAMR