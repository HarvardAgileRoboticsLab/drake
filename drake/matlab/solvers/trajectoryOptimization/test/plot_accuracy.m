close all

options.terrain = RigidBodyFlatTerrain();
options.floating = true;
options.use_bullet = false;

file = fullfile(getDrakePath,'matlab','systems','plants','test','FallingBrickContactPoints.urdf');
plant = RigidBodyManipulator(file,options);
v = plant.constructVisualizer();

q0 = [0
      0
      1.5000
      pi/2
      0
      0];
    
v0 = zeros(6,1);
v0(2) = 5.0;
x0 = [q0;v0];

x0samp = zeros(12,20);
x0samp(:,1) = x0;
for k = 2:20
    x0samp(:,k) = x0 + [.15*randn(6,1); zeros(6,1)];
end

tf=2.0;

ts_plant = TimeSteppingRigidBodyManipulator(plant,0.002,options);

for k = 1:20
sim_traj = ts_plant.simulate([0,tf],x0samp(:,k));
sim_traj = PPTrajectory(foh(sim_traj.getBreaks(), sim_traj.eval(sim_traj.getBreaks())));
true_traj{k} = sim_traj.setOutputFrame(plant.getStateFrame());
end

for k = 1:20
    v.playback(true_traj{k});
end

save('random_true3.mat','true_traj');

% load random_true2.mat

N = 3:2:15;
rms_pos1 = zeros(length(N),20);
rms_pos2 = zeros(length(N),20);

nq=plant.getNumPositions();

for k = 1:20
for i = 1:length(N)
  k
  N(i)
  xtraj1 = variationalBrick(plant,N(i),x0samp(:,k));
  xtraj2 = contactImplicitBrick(plant,N(i),x0samp(:,k));
  %xtraj3 = contactImplicitBrick(plant,N(i),x0,'backward');
  ts1 = xtraj1.getBreaks();
  dt1 = ts1(2)-ts1(1);
  
%   xtraj = contactImplicitBrick(plant,N(i),x0);
%   xtraj_knots = xtraj.eval(xtraj.getBreaks());
%   xtraj2 = PPTrajectory(foh(xtraj.getBreaks()+(0.5*tf/N(i)),xtraj_knots)); % use foh vel
%   xtraj2 = PPTrajectory(foh(xtraj.getBreaks(),xtraj_knots)); % use foh vel

  x1 = xtraj1.eval(ts1);
  x2 = xtraj2.eval(ts1);
  %x3 = xtraj3.eval(ts1);
  q_true = true_traj{k}.eval(ts1);
  %v_true = sim_traj.eval(ts1+.5*dt1);
  rms_pos1(i,k) = rms(vec(x1(1:nq,:)-q_true(1:nq,:)));
  rms_pos2(i,k) = rms(vec(x2(1:nq,:)-q_true(1:nq,:)));
  %rms_pos3(i) = rms(vec(x3(1:nq,:)-q_true(1:nq,:)));
  %rms_vel(i) = rms(vec(x(nq+(1:nq),:)-v_true(nq+(1:nq),:)));
  
  figure(2);
  for j=1:6
    subplot(2,3,j);
    plot(ts1,x1(j,:),'b');
    hold on;
    plot(ts1,x2(j,:),'r');
    %plot(ts1,x3(j,:),'g');
    plot(ts1,q_true(j,:),'k--');
    hold off;
  end
  legend('Variational','Euler','True')
  
end
end

avg_err1 = mean(rms_pos1,2);
avg_err2 = mean(rms_pos2,2);
std1 = std(rms_pos1,1,2);
std2 = std(rms_pos2,1,2);

figure(1)
errorbar(N,avg_err2,std2,'color',[0.8500, 0.3250, 0.0980],'linewidth',2);
hold on;
errorbar(N,avg_err1,std1,'color',[0, 0.4470, 0.7410],'linewidth',2);
legend('Euler','Variational');
xlim([2,16]);

% figure();
% semilogy(N,rms_pos1);
% hold on
% semilogy(N,rms_pos2);

save('variational-rms.mat','N','rms_pos','rms_vel');
