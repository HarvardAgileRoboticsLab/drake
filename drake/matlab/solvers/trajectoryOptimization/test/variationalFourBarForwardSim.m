function [xtraj, kl] = variationalFourBarForwardSim(plant,N,z0,tf) 
addpath('./utils/')

if nargin<1
    options.floating = false;
    options.use_bullet = false;
    file = fullfile(getDrakePath,'examples', 'SimpleFourBar', 'FourBar_JointLimits.urdf');
    plant = RigidBodyManipulator(file,options);    
end
if nargin < 2
    N=21;
end
if nargin<3
    q0 = [4*pi/5; pi/5; 4*pi/5]; %pi/4*ones(3,1);
    v0 = 0*q0; 
    valuecheck(positionConstraints(plant, q0), zeros(6,1), 1e-4);
    x0 = [q0;v0];
end
if nargin<4
    tf = 2;
end

angle_inds = [plant.body(~[plant.body.floating]).position_num]';
angle_inds = angle_inds(angle_inds>0);

h = tf/(N-1);
t = 0:h:tf;

% parameters to pass to initializer
q0 = z0(1:3); v0 = z0(4:6); 
init_params.plant = plant;
init_params.angle_inds = angle_inds;
init_params.h = h;
init_params.q0 = q0;
init_params.v0 = v0;

q = zeros(3, N); q(:,1) = q0; 
kl = zeros(2,N);

% Initial guesses
% kl0 = rand(2,1);
kl0 = [-1.5864; 5.4613];
q1 = q0+0.1*rand(3,1);
zk = [q1; kl0];
[fk, dfk] = initializeDyn(zk, init_params);

% Initialize dyn
while(norm(fk,2) > 1e-6)
    zk = zk - dfk\fk;
    [fk, dfk] = initializeDyn(zk, init_params);
%     disp(fk)
%     disp(zk)
end
fprintf('Cond dF: %s at time %d \r ', cond(dfk), t(2)); 

q(:,2) = zk(1:3);
kl(:,2) = zk(4:5);


%% Simulate forward
params.plant = plant;
params.angle_inds = angle_inds;
params.h = h;
% k = 3;
for k = 3:N
    params.q1 = q(:,k-2);
    params.q2 = q(:,k-1);
    params.kl1 = kl(:,k-1);
    
    zk = [params.q2; params.kl1];
    [fk, dfk] = forwardDyn(zk, params);
    
    while(norm(fk,2) > 1e-6)        
        zk = zk - dfk\fk;
        [fk, dfk] = forwardDyn(zk, params);
%         disp(fk)
%         disp(zk)
        
    end
    fprintf('Cond dF: %s at time %d \r', cond(dfk), t(k));
    q(:,k) = zk(1:3);
    kl(:,k) = zk(4:5);
    
end

%% Validate 

sim_traj = plant.simulate([0. tf], [q0; v0]);
ts = sim_traj.getBreaks(); 
qs = sim_traj.eval(ts); 
% qs_interp = interp1(sim_traj.getBreaks(), qs(1:3,:), t);  


v = zeros(3, N);
for i = 1:N-1
    v(:,i) = qdiff(params, q(:,i), q(:,i+1), h);
end
xtraj = PPTrajectory(foh(t, [q; v])); 

% PLAYBACK
qtraj = PPTrajectory(foh(t, q)); 
vis = plant.constructVisualizer();
qtraj = qtraj.setOutputFrame(vis.getInputFrame); 
vis.playback(qtraj, struct('slider', true))

figure(1); clf; 
for i=1:size(q, 1)  
  subplot(1,size(q,1),i); hold on;
  plot(ts,rad2deg(qs(i,:)),'r');
  plot(t, rad2deg(q(i,:)), 'b--'); 
  hold off;
  legend('MidpointRule', 'ode45')
end


% for i = 1:100
%     zin = 20*rand(size(zi));
%     [f, df] = initializeDyn(zin, init_params);
%
%     df_fd = zeros(size(df));
%     dzin = 1e-6*eye(length(zin));
%     for k = 1:length(zin)
%         df_fd(:,k) = (initializeDyn(zin+dzin(:,k), init_params) - ...
%             initializeDyn(zin-dzin(:,k),init_params))/2e-6;
%     end
%
%     disp('Kinematic Loop Derivative Error:');
%     disp(max(abs(df_fd(:)-df(:))));
% end

% figure(2); clf; 
% for i=1:size(kl, 1)  
%   subplot(1,size(kl,1),i); %hold on;
%   plot(t(2:end), kl(i,(2:end)), 'b--'); 
% %   hold off;
% %   legend('MidpointRule', 'TimeStepping')
% end
