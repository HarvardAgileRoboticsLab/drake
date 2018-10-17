clear; clc; close all;

N = 5;
time = zeros(N,1);

for i = 1:N
    tic;
    [p,xtraj,utraj,ctraj,btraj,psitraj,etatraj,straj,z, ...
    F,info,infeasible_constraint_name,traj_opt] = runStepUp();
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

    data.z = z;
    data.F = F;
    data.info = info;
    data.infeasible_constraint_name = infeasible_constraint_name;
    data.topt = topt;
    save([fname, '.mat'], '-struct', 'data');
    z0 = [];
end
