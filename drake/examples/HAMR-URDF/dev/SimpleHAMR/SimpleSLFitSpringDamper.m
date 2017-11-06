load('sysid_traj.mat')
%% Build Full Single Leg

% options
name = 'FL_scaled';
urdf = fullfile(getDrakePath, 'examples', 'HAMR-URDF', 'dev', 'SimpleHAMR', 'urdf', [name, '.urdf']);
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = false;
options.z_inactive_guess_tol = 0.1;
options.dt = mean(diff(t))/5;

SL = SLTSRBM(urdf, options);
nq = SL.getNumPositions();
nv = SL.getNumVelocities();
nu = SL.getNumInputs();

% Relax joint limits
SL = SL.setJointLimits(-Inf(nq, 1), Inf(nq, 1));
SL = compile(SL); 

%% Build Simple Single Leg

sl_urdf = fullfile(getDrakePath, 'examples', 'HAMR-URDF', 'dev', 'SimpleHAMR', ...
    'urdf','SLSimple_scaled.urdf');

% Initial guess for stiffness and damping 
sfb_inputs = [3; 7]; 
xout = xi(sfb_inputs, :); 
xin = xi(SL.getActuatedJoints(), :); 
Trat = xin'\xout'; 
Trat_sq = Trat*Trat; 

K0 = diag(3*0.05*ones(2,1)) + Trat_sq\diag([1.1556; 1.1556]); 
P0 = diag(0.01*ones(2,1)) + Trat_sq\diag([0.01; 0.01]);

options.k = K0(:);
options.p = P0(:);

SLSimple = SLSimpleRBM(sl_urdf, options);
x0 = SLSimple.getInitialState(); 
nqS = SLSimple.getNumPositions();
nvS = SLSimple.getNumVelocities(); 

%% Fitting

i_start = find(t > 1000, 1, 'first'); 
i_end = i_start + 101; 

params.SL = SL;
params.SLSimple = SLSimple; 
params.t = t(i_start:i_end); 
params.x = xi(1:(nq+nv), i_start:i_end); 
params.u = xi(1+nq+nv+(1:nu), i_start:i_end); 
params.tau = tau(:, i_start:i_end); 
params.sfb_inputs = sfb_inputs; 




