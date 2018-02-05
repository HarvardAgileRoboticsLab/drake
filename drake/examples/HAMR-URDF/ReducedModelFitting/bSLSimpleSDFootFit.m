clear; clc; close all;
addpath('../')
datadir = '~/Dropbox/CurrentWork/FrictionTrajOpt/MatFiles/SysIDFiles/';
dataset = '_40Hz_200V';
data = load([datadir, 'sysid_traj', dataset, '.mat']);

N0 = find(data.t>=500, 1, 'first');
N = 200;
ds = (numel(data.t) - N0)/N; 
IND = N0:ds:numel(data.t) ;      % randomly choose indices
SFB_INDS = [3;7];

NK = 1; 
NP = 1; 

[xS, objval, SL] = SLSimpleIK(data, IND, SFB_INDS); 
fprintf('Average Error: %f mm \r', sqrt(objval/N)); 
[xf,objval,exitflag,infeasible_constraint_name] = SLSimpleFitNLPFoot(data, IND, xS, NK, NP); 

%% Simulate w/ fitted params

sl_urdf = fullfile(getDrakePath, 'examples', 'HAMR-URDF', 'urdf','SLSimple_scaled.urdf');

options.k = reshape(xf(1:4*NK), [], NK);
options.p = reshape(xf(4*NK+(1:4*NP)), [], NP); 

SLSimple = SLSimpleRBM(sl_urdf, options);
v = SLSimple.constructVisualizer();
nqS = SLSimple.getNumPositions();
nvS = SLSimple.getNumVelocities();
x0 = SLSimple.getInitialState();

K = reshape(options.k, nqS, nqS, NK);
P = reshape(options.p, nvS, nvS, NP);

t = data.t; 
tau_traj = PPTrajectory(foh(t, data.tau));
tau_traj = tau_traj.setOutputFrame(SLSimple.getInputFrame());
SLSimple_OL0 = cascade(tau_traj, SLSimple);

xtraj = simulate(SLSimple_OL0, [t(1), t(end)], x0);

%% Validate 
xx = data.x;
xxS = xtraj.eval(t);

xyz_foot = zeros(3, numel(t));
xyz_footS = zeros(3, numel(t));
for i = 1:numel(t)
    
    q = xx(1:SL.getNumPositions(), i);
    qd = xx(SL.getNumPositions()+(1:SL.getNumVelocities()), i);
    kinsol = SL.doKinematics(q, qd);
    xyz_foot(:,i) = SL.forwardKin(kinsol, SL.findLinkId(SL.foot), SL.pf);
    
    qS = xxS(1:nqS, i);
    qdS = xxS(nqS+(1:nvS),i);
    kinsolS = SLSimple.doKinematics(qS, qdS);
    xyz_footS(:,i) = SLSimple.forwardKin(kinsolS, SLSimple.findLinkId('L2'), [0 0 -14.988382167532292]');

end

figure(1); clf; hold on;
titles = {'Leg X (mm)', 'Leg Y (mm)', 'Leg Z (mm)'};
for i = 1:3
    subplot(3,1,i); hold on;
    title(titles{i});
    plot(t, xyz_foot(i,:))
    plot(t, xyz_footS(i,:));
end

%% Sim and save
xtraj_scaled = PPTrajectory(foh(t*1e-3, xxS));
xtraj_scaled = xtraj_scaled.setOutputFrame(xtraj.getOutputFrame());
v.playback(xtraj_scaled, struct('slider', true));

RMS = 2*rms(xyz_foot - xyz_footS, 2)./(var(xyz_foot,[], 2) + var(xyz_footS,[], 2))
save([datadir, 'SpringDamper', dataset '_', num2str(NK), '_', num2str(NP)], 'K', 'P', 'objval'); 

%%
% xd = data.x;
% figure(1); clf; hold on;
% titles = {'Leg X (mm)', 'Leg Y (mm)', 'Leg Z (mm)'};
% for i = 1:3
%     si = subplot(3,1,i); hold on;
%     set(si, 'FontSize', 16); 
%     plot(t*1e-3, , 'LineWidth', 1.5)
%     plot(t*1e-3, , 'LineWidth', 1.5);
%     ylabel(titles{i}, 'FontSIze', 18)
%     lh = legend('Full Model', 'Simple Model')    ;
%     set(lh, 'box', 'off')
%     xlabel('Time (s)')    
%     ylim([-30, 30])
% end
% 
% % rad2deg(rms(bsxfun(@minus, xx, x0),2))
% % rad2deg(rms(xd(ind_plot,:),2))


