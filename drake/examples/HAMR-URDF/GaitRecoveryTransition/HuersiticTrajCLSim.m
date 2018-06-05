function [tt_sol, xx_sol, vv_sol, err_sol] = HuersiticTrajCLSim(hamrWact, hamr, actuators, params)

% unpacks
NCYC = params.NCYC;
RAMPCYC = params.RAMPCYC;
LIFTAMP = params.LIFTAMP;
SWINGAMP = params.SWINGAMP;
gait = params.GAIT;
NPTS = params.NPTS;
X0 = params.X0;
FREQ = params.FREQ;
DC = params.DC_SWING;
DL = params.DC_LIFT;
TYPE = params.TYPE;

% dimensions
nq = hamr.getNumPositions();
nv = hamr.getNumVelocities();
nu = hamr.getNumInputs();
qa = hamr.getActuatedJoints();

tsim = NCYC/FREQ;
t = 0:hamr.timestep:tsim;


%% Load Feed forward inputs
load_dir = '~/Dropbox/GaitRecoveryandTransition/FeedforwardInputs';
ff_inputs = load([load_dir, '/Type', num2str(TYPE) '/', 'FFInputs_Type', num2str(TYPE), '_', num2str(1e3*FREQ), 'Hz.mat']);
ff_inputs = ff_inputs.ffinputs;
DLmat = cellfun(@(x) x.params.DL, ff_inputs);
DCmat = cellfun(@(x) x.params.DC, ff_inputs);
[rind, cind] = find((DLmat == DL).*(DCmat==DC));            % find right input 
vv = repmat(ff_inputs{rind, cind}.vvopt, 1, NCYC);          % copy for NCYC
vtraj = PPTrajectory(foh(mean(diff(ff_inputs{rind, cind}.ttopt))*(0:size(vv, 2)-1), vv));       % 
vtraj = setOutputFrame(vtraj, hamrWact.getInputFrame());
%% Generate Hueristic Trajectories in body frame

% generate trajectory
traj = GenerateHueristicActTraj(gait, FREQ, DC, DL, NPTS, TYPE);
tcyc = traj(:, 1);
qcyc = traj(:,2:end);

% rescale
qcyc(:, 1:2:end) = SWINGAMP*qcyc(:, 1:2:end);
qcyc(:, 2:2:end) = LIFTAMP*qcyc(:, 2:2:end);

% Finite diff to find v_cyc
hhd = (1/FREQ)/NPTS;
e = 1./(2*hhd*ones(NPTS,1));
D1 = spdiags([-e, 0*e, e], [-1, 0, 1], NPTS, NPTS);
D1(1,1) = -3/2/hhd; D1(1,2) = 2/hhd; D1(1,3) =-1/2/hhd;
D1(NPTS,NPTS) = 3/2/hhd; D1(NPTS,NPTS-1) = -2/hhd; D1(NPTS,NPTS-2) = 1/2/hhd;

vcyc = D1*qcyc;
tlabel = {'FL', 'RL', 'FR', 'RR'};
figure(1); clf;
for j =1:nu
    subplot(nu/2, 2, j); hold on;
    plotyy(tcyc, qcyc(:,j), tcyc, vcyc(:,j))
    if rem(j,2) == 0
        title(['L', tlabel{j/2}]);
    else
        title(['S', tlabel{floor(j/2)+1}]);
    end
end

% ramp
tramp = RAMPCYC/FREQ;
ramp = t/tramp; ramp(t >= tramp) = 1;

% create trajectory
ttd = 0:hhd:(tsim-hhd);
xxd = zeros(nq+nv, NPTS*NCYC);
xxa = bsxfun(@times, interp1(t, ramp, ttd), repmat([qcyc, vcyc], NCYC, 1)');
xxd([qa; nq+qa], :) = xxa;
xtrajd = PPTrajectory(foh(ttd, xxd));

%% Simulate Closed Loop

% mirror what's used on the robot
QposL = 0.7*(1/LIFTAMP)^2;
QposS = 7*(1/SWINGAMP)^2;
Qvel = 0.07*(2/(2*pi*FREQ)/(LIFTAMP+SWINGAMP))^2;       % freq normalization means I don't have to rescale Qvel
rho = 0.5*(2/actuators.dummy_bender(1).dp.Vb)^2;

tracking_opt.ctype = 'actlqr';
tracking_opt.rho = rho;
tracking_opt.QposS = QposS;
tracking_opt.QposL = QposL;
tracking_opt.Qvel = Qvel;

LegTracking = HAMRLegTracking(hamrWact, hamr, actuators, vtraj, xtrajd, tracking_opt);

% mimo outputs
output_select(1).system = 1;
output_select(1).output = hamrWact.getOutputFrame.getFrameByName('HamrPosition');
output_select(2).system = 1;
output_select(2).output = hamrWact.getOutputFrame.getFrameByName('HamrVelocity');
output_select(3).system = 2;
output_select(3).output = LegTracking.getOutputFrame();

% simulate
hamr_CL = mimoFeedback(hamrWact, LegTracking, [], [], [], output_select);
xtraj_sim = simulate(hamr_CL, [0, tsim], X0);

tt_sol = xtraj_sim.getBreaks();
yy_sol = xtraj_sim.eval(tt_sol);
xx_sol = yy_sol(1:nq+nv, :);
vv_sol = yy_sol(nq+nv+(1:nu),:);

% calculate tracking error
xx_sol_act = xx_sol(qa,:);
err_sol = xxd(qa,:) - interp1(tt_sol', xx_sol_act', ttd')';

