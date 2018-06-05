clear; clc; close all;
warning('off', 'MATLAB:nargchk:deprecated');
addpath('../', '../urdf/')

%% Parameters of sweep 

freq = linspace(10e-3, 50e-3, 5);           %linspace()*1e-3;
DC = linspace(50, 80, 4);                   % duty cycle for swing
DL = linspace(-25, 75, 5);                  % percent "push" into the ground

Nf = numel(freq);
Ndc = numel(DC);
Ndl = numel(DL); 

NUM = 1; 
MU = 0.51;                                  % this needs to be manually set in @RigidBodyManipulator/ContactConstraints

%% Load Rigid Body Model

urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', 'HAMR_scaledV2.urdf');
save_dir = ['./data_' datestr(date), '/'];
if exist(save_dir, 'dir') == 0
    mkdir(save_dir);
end

% options to change
gait = 'TROT';
NCYC = 25;
RAMPCYC = 10;             % number of cycles to ramp
NPTS = 100;              % number of pts/cycle in desired traj
LIFTAMP = 0.15;          % lift actuator motion (mm)
SWINGAMP = 0.175;        % swing actuator motion (mm)
x0 = zeros(76, 1); x0(3) = 12.69;
dt = 0.4;
TYPE = 1; 

% general rbm options
options.floating = true;
options.collision = true;
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.z_inactive_guess_tol = 0.1;
options.use_bullet = false;
options.terrain = RigidBodyFlatTerrain();

% hamr options
options.stiffness_mult = [1, 1, 1, 1]';     %(swing act, lift act, swing flex, lift flex)

% Build robot + visualizer
hamr = HamrTSRBM(HamrRBM(urdf, options), dt, options);
nq = hamr.getNumPositions();
nv = hamr.getNumVelocities();
nu = hamr.getNumInputs();
qa = hamr.getActuatedJoints();

v = hamr.constructVisualizer();

%% Build Actuators
dp.Vb = 225;
dp.Vg = 0;

nact = 8;
hr_actuators = HamrActuators(nact, {'FLsact', 'FLlact', 'RLsact', 'RLlact', ...
    'FRsact', 'FRlact', 'RRsact', 'RRlact'}, [], dp);

%% Connect system

%connections from actuators to hamr
hr_actuators = hr_actuators.setOutputFrame(hamr.getInputFrame());
connection1(1).from_output = hr_actuators.getOutputFrame();
connection1(1).to_input = hamr.getInputFrame();

% connections from hamr to actuators
hamr_out = hamr.getOutputFrame();
act_in = hr_actuators.getInputFrame();
act_in = act_in.replaceFrameNum(2, hamr_out.getFrameByName('ActuatorDeflection'));
hr_actuators = hr_actuators.setInputFrame(act_in);

connection2(1).from_output = hamr_out.getFrameByName('ActuatorDeflection');
connection2(1).to_input = act_in.getFrameByName('ActuatorDeflection');

% mimo inputs
input_select(1).system = 1;
input_select(1).input = act_in.getFrameByName('DriveVoltage');

% mimo outputs
output_select(1).system = 2;
output_select(1).output = hamr_out.getFrameByName('HamrPosition');
output_select(2).system = 2;
output_select(2).output = hamr_out.getFrameByName('HamrVelocity');

hamrWact = mimoFeedback(hr_actuators, hamr, connection1, connection2, ...
    input_select, output_select);

%% Simulate for all freq, dc
params.NCYC = NCYC;
params.RAMPCYC = RAMPCYC; 
params.LIFTAMP = LIFTAMP;
params.SWINGAMP = SWINGAMP;
params.GAIT = gait;
params.NPTS = NPTS;
params.TYPE = TYPE;
params.STIFFNESS_MULT = options.stiffness_mult; 
params.MU = MU; 

params.X0 = x0;

cl_sols = cell(Ndc, Ndl);

for i = 1:Nf

    for j = 1:Ndc
        for k = 1:Ndl
            params.FREQ = freq(i);
            params.DC_SWING = DC(j);
            params.DC_LIFT = DL(k);
            sol_struct.params = params;
            try
                fprintf('Simulating %d Hz at %d DC Swing and %d DC Lift \r', freq(i)*1000, DC(j), DL(k))
                tic;
                [tt_sol, xx_sol, vv_sol, err_sol] = HuersiticTrajCLSim(hamrWact, hamr, hr_actuators, params);
                fprintf('This took %f s \r', toc)
                sol_struct.tt_sol = tt_sol;
                sol_struct.xx_sol = xx_sol;
                sol_struct.vv_sol = vv_sol;
                sol_struct.err_sol = err_sol;
                                
                %% calculate leg position
                pfCL = zeros([numel(tt_sol), size(lp_b')]);
                vfCL = zeros([numel(tt_sol), size(lp_b')]);
                legs = {'FLL4', 'RLL4', 'FRL4', 'RRL4'};
                
                for l = 1:numel(tt_sol)
                    qCL = xx_sol(1:nq, l);
                    qdCL = xx_sol(nq+1:2*nq, l);
                    kinsolCL = hamr.doKinematics(qCL, qdCL, struct('compute_gradient', true));
                    for m = 1:size(lp_b,1)
                        [pfCL(l,:,m), J] = hamr.forwardKin(kinsolCL, hamr.findLinkId(legs{m}), lp_b(m,:)'); %,fkopt);
                        vfCL(l,:,m) = J*qdCL;
                    end
                end
                
                sol_struct.pfCL = pfCL;
                sol_struct.vfCL = vfCL;
                
                
            catch ME
                disp('')
                if strcmp(ME.identifier, 'Simulink:SFunctions:SFcnErrorStatus')
                    sol_struct.tt_sol = [];
                    sol_struct.xx_sol = [];
                    sol_struct.vv_sol = [];
                    sol_struct.err_sol = [];
                    sol_struct.pfCL = [];
                    sol_struct.vfCL = [];
                    
                else
                    throw(ME)
                end
            end            
            cl_sols{j, k} = sol_struct;            
        end
    end
    disp('Saving...')
    save([save_dir, 'ModelSim_Type', num2str(TYPE), '_'  gait, '_', num2str(1e3*freq(i)), 'Hz_', num2str(NUM) ], ...
        'cl_sols')
    
end



