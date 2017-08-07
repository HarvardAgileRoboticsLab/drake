clear; clc; close all;
load Trot1Hz.mat
load Trot_1Hz_175V.mat
te = Trot1Hz.Time + 0.228;
t0 = xtraj_scaled.getBreaks();
x0 = xtraj_scaled.eval(t0);
tm = t0(t0>3);
xm = x0(:, t0>3);
% tm = tm - tm(1);
RPY = [Trot1Hz.RollDeg', Trot1Hz.PitchDeg', Trot1Hz.YawDeg'];  
for i = 1:3
    xyze(i,:) = Trot1Hz.XYZCOM(i,:)*1e3 - polyval(polyfit(te,Trot1Hz.XYZCOM(i,:)*1e3, 0), te);
    rpye(:,i) = RPY(:,i) - polyval(polyfit(te',RPY(:,i), 0), te');
end
% xyze = Trot1Hz.XYZCOM*1e3;
xyzm = xm(1:3, :);
rpym = xm(4:6, :); 

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
titles = {'Roll', 'Pitch', 'Yaw'};
for i = 1:size(xyze,1)
    si = subplot(3,1,i); hold on; title(titles{i}, 'FontSize', 20);
    set(si, 'fontsize', 14);
    plot(tm, rad2deg(rpym(i,:) - mean(rpym(i,:))), te, rpye(:,i))
    xlim([2, 8]);
    if i == 2
        ylabel('Orientation(deg)', 'FontSize', 16)
    end
    if i == 3
        xlabel('Time(s)', 'FontSize', 16);
    end
end
lh = legend('Model', 'Experimental');
set(lh, 'box', 'off')

figure(3); clf; hold on;
lege = [Trot1Hz.FLfootXYZ(3,:); Trot1Hz.RLfootXYZ(3,:); Trot1Hz.FRfootXYZ(3,:); Trot1Hz.RRfootXYZ(3,:)]*1e3;
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
