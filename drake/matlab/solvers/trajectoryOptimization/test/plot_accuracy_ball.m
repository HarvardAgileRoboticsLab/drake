close all

options.terrain = RigidBodyFlatTerrain();
options.floating = true;
options.use_bullet = true;

file = fullfile(getDrakePath,'matlab','systems','plants','test','ball.urdf');
plant = RigidBodyManipulator(file,options);
v = plant.constructVisualizer();

q0 = [0
      0
      0.5000
      0
      0
      0];
    
v0 = zeros(6,1);
v0(2) = 3.0;
v0(4) = 10.0;
x0 = [q0;v0];

tf=2.0;

ts_plant = TimeSteppingRigidBodyManipulator(plant,0.0005,options);

sim_traj = ts_plant.simulate([0,tf],x0);
sim_traj = PPTrajectory(foh(sim_traj.getBreaks(), sim_traj.eval(sim_traj.getBreaks())));
true_traj = sim_traj.setOutputFrame(plant.getStateFrame());
v.playback(true_traj);

N = 5:2:15;
rms_pos1 = zeros(length(N),1);
rms_pos2 = zeros(length(N),1);

nq=plant.getNumPositions();

for i = 1:length(N)
  N(i)
  xtraj1 = variationalBrick(plant,N(i),x0);
  xtraj2 = contactImplicitBrick(plant,N(i),x0);
  ts1 = xtraj1.getBreaks();
  dt1 = ts1(2)-ts1(1);
  
%   xtraj = contactImplicitBrick(plant,N(i),x0);
%   xtraj_knots = xtraj.eval(xtraj.getBreaks());
%   xtraj2 = PPTrajectory(foh(xtraj.getBreaks()+(0.5*tf/N(i)),xtraj_knots)); % use foh vel
%   xtraj2 = PPTrajectory(foh(xtraj.getBreaks(),xtraj_knots)); % use foh vel

  x1 = xtraj1.eval(ts1);
  x2 = xtraj2.eval(ts1);
  q_true = true_traj.eval(ts1);
  rms_pos1(i) = rms(vec(x1(1:nq,:)-q_true(1:nq,:)));
  rms_pos2(i) = rms(vec(x2(1:nq,:)-q_true(1:nq,:)));
  
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

figure(1)
plot(N,rms_pos2,'color',[0.8500, 0.3250, 0.0980],'linewidth',2);
hold on;
plot(N,rms_pos1,'color',[0, 0.4470, 0.7410],'linewidth',2);
legend('Euler','Variational');
%xlim([2,16]);

% figure();
% semilogy(N,rms_pos1);
% hold on
% semilogy(N,rms_pos2);

save('variational-rms.mat','N','rms_pos','rms_vel');
