clear; clc; close all;

load HamrVarViz_175V_1Hz_Trot.mat
ttv = tt; 
Vactv = Vact; 
yyv = yy; 
load HamrViz_175V_1Hz_Trot.mat

figure(1); clf; hold on; 
for i = 1:6
    subplot(2,3,i); hold on; 
    plot(tt, yy(i,:))
    subplot(2,3,i); hold on; 
    plot(ttv, rad2deg(yyv(i,:)), '--')
    legend('Drake Const', 'New Const')
end

figure(2); clf; hold on; 
for i = 1:8
    subplot(2,4,i); hold on; 
    plot(tt, Vact(i,:))
    subplot(2,4,i); hold on; 
    plot(ttv, Vactv(i,:), '--')
    legend('Drake Const', 'New Const')
end

figure(3); clf; hold on; 
for i = 1:8
    subplot(2,4,i); hold on; 
    plot(tt, yy(76+i,:))
    subplot(2,4,i); hold on; 
    plot(ttv, yyv(100+i,:), '--')
    legend('Drake Const', 'New Const')
end

% for i = 7:38
%     figure(i); clf; hold on;
%     plot(tt, yy(i,:))
%     plot(ttv, yyv(i,:), '--')
% end

tilefigs;

