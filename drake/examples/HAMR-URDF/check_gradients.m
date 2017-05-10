clear; clc; close all;

%% BUILD ROBOT
pkg_path = '/home/nddoshi/dev/drake/drake/examples/HAMR-URDF/';
urdf = [pkg_path, 'SL_assem6/urdf/HAMR_scaled.urdf'];

% options
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.z_inactive_guess_tol = 1;
options.dt = 2;  %time step
options.use_bullet = false;

% floating (gnd contact) or in air (not floating)
ISFLOAT = true;
% ISFLOAT = true;

if ISFLOAT
    options.floating = ISFLOAT;
    options.collision = ISFLOAT;
    options.terrain = RigidBodyFlatTerrain();
else
    options.floating = ISFLOAT;
    options.collision = ISFLOAT;
    options.terrain = [];
end

% Build robot + visualizer
hamr = Hamr(urdf, options);
hamr = compile(hamr);

Nx = hamr.getNumStates();
Nu = hamr.getNumInputs();

CHECKTOL = 1e-6;
 
%% For RIGID BODY FORCE ELEMENTS
n = 10;

hmanip = hamr.getManipulator();
forces = hmanip.force;
for j = 1:numel(forces)
    passed = 0;
    for i = 1:n
        q = rand(Nx/2,1);
        qd = rand(Nx/2, 1); 
               
        [~, dfext_dx] = forces{j}.computeSpatialForce(hmanip, q, qd); 
        dfext_dx = reshape(dfext_dx, [], Nx); 
        f = @(q, qd) forces{j}.computeSpatialForce(hmanip, q, qd); 
        
        [~, dfext_dx2] = geval(f, q, qd, struct('grad_method','numerical'));
        passed = passed + valuecheck(dfext_dx,dfext_dx2, CHECKTOL);        
    end
    if sum(passed == n); disp(['gradients for force ' num2str(j) ' are correct']); end
end

%% For RIGID BODY LOOPS
n = 10;

hmanip = hamr.getManipulator();
loops = hmanip.loop;
for j = 1:numel(loops)
    passed = 0;
    for i = 1:n
        q = rand(Nx/2,1);
        qd = rand(Nx/2, 1); 
               
        [~, dfext_dx] = loops(j).computeSpatialForce(hmanip, q, qd); 
        dfext_dx = reshape(dfext_dx, [], Nx); 
        f = @(q, qd) loops(j).computeSpatialForce(hmanip, q, qd); 
        
        [~, dfext_dx2] = geval(f, q, qd, struct('grad_method','numerical'));
        passed = passed + valuecheck(dfext_dx,dfext_dx2, CHECKTOL);        
    end
    if sum(passed == n); disp(['gradients for loop ' num2str(j) ' are correct']); end
end


%% For HAMR OUTPUT

n = 10;
passed = 0;
for j = 1:n
    t = rand(1);
    x = rand(Nx,1);
    u = rand(Nu,1);
    
    [~, dfext_dx] = hamr.output(t,x,u);
    f = @(t, x,u) hamr.output(t,x,u);
    
    [~, dfext_dx2] = geval(f, t, x, u, struct('grad_method','numerical'));
    passed = passed + valuecheck(dfext_dx,dfext_dx2, CHECKTOL);
end

if sum(passed == n); disp('Hamr ouput gradients are correct'); end