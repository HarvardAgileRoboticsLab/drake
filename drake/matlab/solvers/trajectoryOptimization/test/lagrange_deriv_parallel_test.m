clear; clc; close all; 
addpath('./utils/')

urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', 'HAMR_scaled.urdf');

options.ignore_self_collisions = true;
options.use_bullet = false;
options.floating = true;
options.collision = false;


obj = RigidBodyManipulator(urdf, options);
nq = obj.getNumPositions();


N = 100;
for i = 1:N
    x0 = rand(2*nq,1); 
    q0 = x0(1:nq);
    v0 = x0(nq+1:2*nq); 
    
    tic;
    [D1L,D2L,D1D1L,D1D2L,D2D2L,B,dBdq] = LagrangianDerivs(obj,q0,v0);
    t = toc; disp(['Serial : ', num2str(t)])
    
%     tic;
    [D1LP,D2LP,D1D1LP,D1D2LP,D2D2LP,BP,dBdqP] = LagrangianDerivsParallel(obj,q0,v0);
%     t = toc; disp(['Parallel : ', num2str(t)])
    
    valuecheck(D1D1L, D1D1LP, 1e-6)
end
