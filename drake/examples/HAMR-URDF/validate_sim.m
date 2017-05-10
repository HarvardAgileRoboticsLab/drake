% clear; clc; close all;
load Mets.mat
load Trot_1Hz_175V.mat
te = Mets.Time + 0.228;
t0 = xtraj_scaled.getBreaks();
x0 = xtraj_scaled.eval(t0);
tm = t0(t0>3);
xm = x0(:, t0>3);
% tm = tm - tm(1);
for i = 1:3
    xyze(i,:) = Mets.XYZCOM(i,:)*1e3 - polyval(polyfit(te,Mets.XYZCOM(i,:)*1e3, 0), te);
end
% xyze = Mets.XYZCOM*1e3;
xyzm = xm(1:3, :);

figure(1); clf; hold on;
titles = {'Com-X', 'Com-Y', 'Com-Z'};
for i = 1:size(xyze,1)
    si = subplot(3,1,i); hold on; title(titles{i}, 'FontSize', 20);
    set(si, 'fontsize', 14);
    plot(tm, xyzm(i,:) - mean(xyzm(i,:)), te, xyze(i,:))
    xlim([2, 8]);
    if i == 2
        ylabel('Position(mm)', 'FontSize', 16)
    end
    if i == 3
        xlabel('Time(s)', 'FontSize', 16);
    end
end
lh = legend('Model', 'Experimental');
set(lh, 'box', 'off')


figure(2); clf; hold on;
lege = [Mets.FLfootXYZ(3,:); Mets.RLfootXYZ(3,:); Mets.FRfootXYZ(3,:); Mets.RRfootXYZ(3,:)]*1e3;
titles = {'FL Leg', 'RL Leg', 'FR Leg', 'RR Leg'};
for i = 1:size(lege,1)
    si = subplot(2,2,i); hold on; title(titles{i}, 'fontsize', 20);
    set(si, 'fontsize', 14);
    plot(tm, phi(i, t0>3), te, lege(i,:))
    xlim([2, 8]);
    ylim([0, 6]); 
    if i == 3 || i == 4
        xlabel('Time(s)', 'FontSize', 16);
    end
    if i == 1 || i == 3
        ylabel('Height(mm)', 'fontsize', 16)
    end
end
lh = legend('Model', 'Experimental');
set(lh, 'box', 'off')
