clear; clc; close all;

%% Build (open-loop) control input

fname = 'TrajOpt_26-Oct-2017_6';
% traj_COM = load(fname);
traj_full = load([fname, '_fullRobot']);

% Build transmission trajectory
xtrajTrans = traj_full.xtraj();
tt = xtrajTrans.FL_scaled.getBreaks();
hh = mean(diff(tt));

xxFL = xtrajTrans.FL_scaled.eval(tt);
xxRL = xtrajTrans.RL_scaled.eval(tt);
xxFR = xtrajTrans.FR_scaled.eval(tt);
xxRR = xtrajTrans.RR_scaled.eval(tt);

nqt = size(xxFL,1)/2;
nut = 2;        % inputs per transmission 
nct = 1;        % contacts per transmission 
ndt = 4;        % pyramidal friction 
nlt = 18;       % loop const per trans

% Full trajectory
xx = [xxFL(1:nqt, :); xxRL(1:nqt, :); xxFR(1:nqt, :); xxRR(1:nqt, :);
   xxFL(nqt+(1:nqt), :); xxRL(nqt+(1:nqt), :); xxFR(nqt+(1:nqt), :); xxRR(nqt+(1:nqt), :)];

% Build input trajectory (included contact and const. forces)
utrajTrans = traj_full.utraj;
uuFL = utrajTrans.FL_scaled.eval(tt+hh/2);
uuRL = utrajTrans.RL_scaled.eval(tt+hh/2);
uuFR = utrajTrans.FR_scaled.eval(tt+hh/2);
uuRR = utrajTrans.RR_scaled.eval(tt+hh/2);

uu = [uuFL(1:nut,:); uuRL(1:nut,:); uuFR(1:nut,:); uuRR(1:nut,:);];
ll = [uuFL(nut+(1:nlt),:); uuRL(nut+(1:nlt),:); uuFR(nut+(1:nlt),:); uuRR(nut+(1:nlt),:);];
cc = [uuFL(nut+nlt+(1:nct),:); uuRL(nut+nlt+(1:nct),:); uuFR(nut+nlt+(1:nct),:); uuRR(nut+nlt+(1:nct),:)];
bb = [uuFL(nut+nlt+nct+(1:nct*ndt),:); uuRL(nut+nlt+nct+(1:nct*ndt),:); uuFR(nut+nlt+nct+(1:nct*ndt),:); uuRR(nut+nlt+nct+(1:nct*ndt),:)];

%% Build Model
options.terrain = RigidBodyFlatTerrain();
options.ignore_self_collisions = true;
options.collision_meshes = false;
options.use_bullet = false;
options.floating = false;
options.collision = true;
% options.dt = 1; 

hamr = HamrRBM(options);
v = hamr.constructVisualizer();

%% Simulate

tsim = tt(end)/8;
utraj = PPTrajectory(foh(tt, [uu; ll; cc; bb]));
utraj = utraj.setOutputFrame(hamr.getInputFrame);
xtraj = PPTrajectory(foh(tt, xx));
hamr_OL = cascade(utraj, hamr);

x0 = xx(:,1);
[tf, err_str] = valuecheck(positionConstraints(hamr,double(x0)),zeros(72,1),1e-6);

if tf
    disp('Valid initial condition: simulating...')
    tic;
    xtraj_sim = simulate(hamr_OL, [0 tsim], double(x0));
    tlcp = toc;
    xtraj_scaled = PPTrajectory(foh(xtraj_sim.getBreaks()*1e-3, xtraj_sim.eval(xtraj_sim.getBreaks())));
    xtraj_scaled = xtraj_scaled.setOutputFrame(xtraj_sim.getOutputFrame());
    fprintf('It took %fs to simulate %fs of realtime. \nThats %fx \n', ...
        tlcp, tsim/1000, 1000*tlcp/tsim)
    options.slider = true;
    %     xtraj.tt = xtraj.tt/1000;
    v.playback(xtraj_scaled, options);
else
    disp('invalid initial condition...')
end
