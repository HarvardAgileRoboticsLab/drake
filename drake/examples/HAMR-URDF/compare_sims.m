clc; clear; close all; 
load Trot_1Hz_175V.mat
ttOG = xtraj_scaled.getBreaks();
xxOG = xtraj_scaled.eval(xtraj_scaled.getBreaks());
phiOG = phi;

load Trot_1Hz_175V_2.mat
tt = xtraj_scaled.getBreaks();
xx = xtraj_scaled.eval(xtraj_scaled.getBreaks());

figure(1); clf; hold on;
subplot(3,1,1); hold on; 
plot(ttOG, xxOG(1,:),tt, xx(1,:))
subplot(3,1,2); hold on; 
plot(ttOG, xxOG(2,:),tt, xx(2,:))
subplot(3,1,3); hold on; 
plot(ttOG, xxOG(3,:),tt, xx(3,:))

figure(2); clf; hold on;
subplot(2,2,1); hold on; 
plot(ttOG, phiOG(1,:),tt, phi(1,:))
subplot(2,2,2); hold on; 
plot(ttOG, phiOG(2,:),tt, phi(2,:))
subplot(2,2,3); hold on; 
plot(ttOG, phiOG(3,:),tt, phi(3,:))
subplot(2,2,4); hold on; 
plot(ttOG, phiOG(4,:),tt, phi(4,:))
% 
% LIA = ismembertol(xx, xxOG, 1e-4);
% sum(sum(LIA))/numel(LIA)
% spy(LIA)