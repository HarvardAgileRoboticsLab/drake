function [] = GenerateActuatorDeflections()


% Leg Positions
Nz = 5;
zz = linspace(-1.125, 0.375, Nz);

% Leg Deltas
Nd = 3;
dz = linspace(0.1, 0.3, Nd);

delta_lift_mat = cell(Nz, Nd); 

% NPTS
NPTS = 10;

% SAVE OPTIONS
SAVE = 0; 

urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', 'HAMR_scaledV2.urdf');

% options
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.z_inactive_guess_tol = 0.1;
options.use_bullet = false;

% options to change
dt = 0.2;
options.dt = dt;
options.floating = false;
options.collision = false;

% Build robot + visualizer
hamr = HamrRBM(urdf, options);


figure(1); clf; 
for j = 1:Nd
    for i = 1:Nz
        delta_lift = compute_actuator_deflection(hamr, zz(i), dz(j), NPTS);
        if SAVE; save_data(delta_lift, zz(i), dz(j)); end
        % Add to matrix
        delta_lift_mat{i,j} = delta_lift;
        
        % Plot
        subplot(Nd, Nz, (j-1)*Nz + i); hold on; 
        title(sprintf('LegZ: %f, DeltaZ: %f', zz(i), dz(j))) 
        plot(zz(i) + linspace(0, dz(j), NPTS), abs(delta_lift_mat{i, j}), '*-')         
        xlabel('Leg Height (mm)'); 
        ylabel('Lift Actuator Position (mm)')
        xlim([-1.125, .5])
        ylim([0, 0.12])
%         axis equal;
    end
end
end



function dlift = compute_actuator_deflection(hamr, z0, dz, NPTS)


nq = hamr.getNumPositions();
all_act = hamr.getActuatedJoints();
lift_act = all_act(2:2:end);

q0 = zeros(nq, 1);
dz = linspace(0, dz, NPTS);
dlift = zeros(NPTS, numel(lift_act));

v = hamr.constructVisualizer(); 

for i = 1:NPTS
    qi = HAMR_IK(hamr, 0, z0+dz(i), q0);
    dlift(i,:) = qi(lift_act);
    q0 = qi; 
%     v.inspector([qi; 0*qi])
end

end

function [] = save_data(dlift, legz, dz)

    save_dir = '/home/nddoshi/Dropbox/GaitRecoveryandTransition/data/LegImpedanceData/ActuatorDeflections/';
    fname = sprintf('legz_%g_dz_%g.mat', legz, dz);
    FL = dlift(:,1); RL = dlift(:,2); FR = dlift(:,3); RR = dlift(:,4);
    save([save_dir, fname], 'FL', 'RL', 'FR', 'RR', 'legz', 'dz')

end
