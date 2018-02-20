function plotting(t0, tf, tf_phi, xtraj, xtraj_discrete, xtraj_tvlqr, xtraj_phi)
    t = linspace(t0, tf, 100);

    % define more informative time for phase-based controller
    tphi = linspace(t0, tf_phi, 100);

    % evaluate trajectories
    xtraj_discrete_eval = xtraj_discrete.eval(t);
    xtraj_tvlqr_eval = xtraj_tvlqr.eval(t);
    xtraj_phi_eval = xtraj_phi.eval(tphi);
    xtraj_eval = xtraj.eval(t);

    % define states for labels
    states = {'x', 'y', 'z', 'q1', 'q2', 'q3', 'q4', 'vx', 'vy', 'vz', 'w1', 'w2', 'w3'};

    % plot all states for all controller trajectories
    figure('Position', [100 100 1000 800]);
    for i = 1:12
        subplot(4, 4, i)
        plot(t, xtraj_eval(i,:), 'r'), hold on
        plot(t, xtraj_discrete_eval(i,:), 'k--'),
        plot(t, xtraj_tvlqr_eval(i,:), 'k'),
        plot(tphi, xtraj_phi_eval(i,:), 'b'), hold off
        title(states{i}, 'Fontsize', 14);
        set(gca, 'Xlim', [t0 tf_phi]);
        ytickformat('%.2f');
    end
    subplot(4, 4, [14 15])
    plot(t, xtraj_eval(13,:), 'r'), hold on
    plot(t, xtraj_discrete_eval(13,:), 'k--'),
    plot(t, xtraj_tvlqr_eval(13,:), 'k'),
    plot(tphi, xtraj_phi_eval(13,:), 'b'), hold off
    title(states{13}, 'Fontsize', 14);
    set(gca, 'Xlim', [t0 tf_phi]);
    legend('nominal', 'discrete TVLQR', 'TVLQR', 'phi-TVLQR', 'Location', 'bestoutside');
    ytickformat('%.2f');

    % plot 3D positions for all controller trajectories
    figure
    plot3(xtraj_eval(1,:), xtraj_eval(2,:), xtraj_eval(3,:), 'r'); hold on;
    plot3(xtraj_discrete_eval(1,:), xtraj_discrete_eval(2,:), xtraj_discrete_eval(3,:), 'k--'),
    plot3(xtraj_tvlqr_eval(1,:), xtraj_tvlqr_eval(2,:), xtraj_tvlqr_eval(3,:), 'k'),
    plot3(xtraj_phi_eval(1,:), xtraj_phi_eval(2,:), xtraj_phi_eval(3,:), 'b'), hold off
    %set(gca, 'xlim', [-6 6], 'ylim', [-6 6], 'zlim', [-6 6]);
    title('Plane trajectory in space');
    xlabel('x'), ylabel('y'), zlabel('z');
    legend('nominal', 'discrete TVLQR', 'TVLQR', 'phi-TVLQR', 'Location', 'northeast');

    % plot distances from the collision point for all trajectories
    figure
    nom_dist = sqrt((xtraj_eval(1, :) - xtraj_eval(1, end)).^2 + ...
                    (xtraj_eval(2, :) - xtraj_eval(2, end)).^2 + ...
                    (xtraj_eval(3, :) - xtraj_eval(3, end)).^2);
    tvlqr_dist = sqrt((xtraj_tvlqr_eval(1, :) - xtraj_eval(1, end)).^2 + ...
                    (xtraj_tvlqr_eval(2, :) - xtraj_eval(2, end)).^2 + ...
                    (xtraj_tvlqr_eval(3, :) - xtraj_eval(3, end)).^2);
    phi_dist = sqrt((xtraj_phi_eval(1, :) - xtraj_eval(1, end)).^2 + ...
                    (xtraj_phi_eval(2, :) - xtraj_eval(2, end)).^2 + ...
                    (xtraj_phi_eval(3, :) - xtraj_eval(3, end)).^2);
    discrete_dist = sqrt((xtraj_discrete_eval(1, :) - xtraj_eval(1, end)).^2 + ...
                    (xtraj_discrete_eval(2, :) - xtraj_eval(2, end)).^2 + ...
                    (xtraj_discrete_eval(3, :) - xtraj_eval(3, end)).^2);
    plot(t, tvlqr_dist, tphi, phi_dist, t, discrete_dist)
    legend('tvlqr', 'phi-tvlqr', 'discrete tvlqr')
end