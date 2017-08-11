clear; clc; close all;


errors = dir('*.mat');

dt_vec = zeros(numel(errors,1));
rms_err_vec = zeros(numel(errors, 1)); 
for i = 1:numel(errors)
    load(errors(i).name)    
    dt_vec(i) = dt;
    rms_err_vec(i) = mean(rms_err(1:3)); 
end
pft = polyfit(log(dt_vec), log(rms_err_vec), 1)

figure(1); clf;
set(gca, 'xdir', 'reverse')
plot(log(dt_vec), log(rms_err_vec), 'r*'); hold on;
plot(log(dt_vec), polyval(pft, log(dt_vec))); 
% loglog(dt_vec, 2*dt_vec - exp(-1.2818), 'k'); hold on;

% loglog(dt_vec, polyval(polyfit(dt_vec, rms_err_vec, 1), dt_vec)); 