clear; clc; close all;

N = 5;
time = zeros(N,1);

for i = 1:N
    tic;
    [p,xtraj,utraj,ctraj,btraj,psitraj,etatraj,straj,z, ...
    F,info,infeasible_constraint_name,traj_opt] = variationalTrajOpt();
    topt = toc;
    time(i) = topt;
    disp(time(i));
    fname = ['SFTimedTrial_', num2str(i)];
    
    % final state
    data.t = xtraj.getBreaks();
    data.x = xtraj.eval(data.t);
    data.u = utraj.eval(data.t);
    data.c = ctraj.eval(data.t);
    data.b = btraj.eval(data.t);
    data.psi = psitraj.eval(data.t);
    data.eta = etatraj.eval(data.t);
    data.s = straj.eval(data.t);
    
    % initial state
%     data.t0 = z0(traj_opt.h_inds);
%     data.x0 = z0(traj_opt.x_inds);
%     data.u0 = z0(traj_opt.u_inds);
%     data.c0 = z0(traj_opt.c_inds);
%     data.b0 = z0(traj_opt.b_inds);
%     data.psi0 = z0(traj_opt.psi_inds);
%     data.eta0 = z0(traj_opt.eta_inds);
%     data.s0 = z0(traj_opt.s_inds);       
    
    data.z = z;
    data.F = F;
    data.info = info;
    data.infeasible_constraint_name = infeasible_constraint_name;
    data.topt = topt;
    save([fname, '.mat'], '-struct', 'data');
end

% v = p.constructVisualizer(); 
% v.playback(xtraj, struct('slider', true))
