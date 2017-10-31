clear; clc; close all; 
% Drive parameters
Vdrive = 225;   % drive voltage (V)
f = 1;          % drive frequency (Hz)

% Actuator Thickness
t_act = 320e-6;     % m

% time vector
dt = 0.005;             %s
t = 0:dt:(1/f); 
V = 0.5*(Vdrive*sin(2*pi*f*t) + Vdrive); 

% actuator
actuator =  PZTBender2(); 
% actuator = zhi_actuator(t_act); 

% Dynamic model 
delta = zeros(size(t));     % actuator deflection (m)
F_top = zeros(size(t));     % top plate force (N)
F_bot = zeros(size(t));     % botton plate force (N)
E_ave = zeros(size(t));     % botton plate force (N)
K = zeros(size(t)); 

for i = 1:numel(t)
    [F_top(i+1), F_bot(i+1), K(i+1), E_ave(i+1)] = actuator.output(0, delta(i), V(i));
    delta(i+1) = (F_top(i+1) - F_bot(i+1))/K(i+1); 
end

delta = delta(2:end); 
F_top = F_top(2:end);
F_bot = F_bot(2:end);
K = K(2:end);
E_ave = E_ave(2:end); 


%% 

f1 = figure(1); clf; 
set(f1, 'color', 'w')

s1 = subplot(1,3,1); hold on; 
set(s1, 'fontsize', 20)
plot(t, (F_top - F_bot)*1e3, 'k--', 'LineWidth', 3);
% plot(t, Fb*ones(numel(t),1)*1e3, 'k' ,'LineWidth', 3);
% plot(t, -Fb*ones(numel(t),1)*1e3, 'k', 'LineWidth', 3)
% plot(t', Fb_simple*ones(numel(t),1)*1e3, 'b','LineWidth', 3); 
% plot(t', -Fb_simple*ones(numel(t),1)*1e3, 'b', 'LineWidth', 3)
ylabel('Actuator Force (mN)', 'FontSize', 22)

s2 = subplot(1,3,2); hold on;
set(s2, 'fontsize', 20)
plot(t, delta*1e6, 'k--', 'LineWidth', 3);
% % plot(t, 0.5*df*ones(numel(t),1)*1e6, 'k','LineWidth', 3)
% plot(t, -0.5*df*ones(numel(t),1)*1e6, 'k', 'LineWidth', 3)
% plot(t, 0.5*df_simple*ones(numel(t),1)*1e6, 'b', 'LineWidth', 3); 
% plot(t, -0.5*df_simple*ones(numel(t),1)*1e6, 'b', 'LineWidth', 3)
ylabel('Deflection (\mu m)', 'FontSize', 22)
xlabel('Time(s)', 'FontSize', 22)

s3 = subplot(1,3,3); hold on;
set(s3, 'FontSize', 20)
p1 = plot(t, K, 'k--', 'LineWidth', 3);
% p2 = plot(t, 2*Fb/df*ones(numel(t),1), 'k', 'LineWidth', 3);
% p3 = plot(t, 2*Fb_simple/df_simple*ones(numel(t),1), 'b', 'LineWidth', 3); 
ylabel('Stiffness (N/m)', 'FontSize', 22)
% l1 = legend([p1, p2, p3], 'Time Varying Full Model', ...
%     'Full Model', 'Simple Model');
% set(l1, 'box', 'off')

figure; 
plot(delta, E_ave)

