% clear; clc; close all;
% 
% N = 30;
% time = zeros(N,1);
% 
% for i = 21:N
%     tic;
%     [p,xtraj,utraj,ctraj,btraj,psitraj,etatraj,straj,z, ...
%     F,info,infeasible_constraint_name,traj_opt] = runStepUp();
%     topt = toc;
%     time(i) = topt;
%     disp(time(i));
%     fname = ['LDTimedTrial_', num2str(i)];
% 
%     % final state
%     data.t = xtraj.getBreaks();
%     data.x = xtraj.eval(data.t);
%     data.u = utraj.eval(data.t);
%     data.c = ctraj.eval(data.t);
%     data.b = btraj.eval(data.t);
%     data.psi = psitraj.eval(data.t);
%     data.eta = etatraj.eval(data.t);
%     data.s = straj.eval(data.t);
% 
%     data.z = z;
%     data.F = F;
%     data.info = info;
%     data.infeasible_constraint_name = infeasible_constraint_name;
%     data.topt = topt;
%     save([fname, '.mat'], '-struct', 'data');
%     z0 = [];
% end
%% Plot Results
clear; close all;
all_data = dir('./TimingResults/*.mat');

N = numel(all_data);

for i = 1:N
    
    datai = load(['./TimingResults/', all_data(i).name]);
    if isempty(datai.infeasible_constraint_name)
        disp([all_data(i).name, ' converged']);
        t(i) = datai.topt;
        s(i,:) = datai.sXSSLA
        s_inf(i) = norm(datai.s, inf);
        s_one(i) = norm(datai.s, 1)/numel(s);
        s_freq(i) = sum((datai.s >=1e-4))/numel(s);
    else
        disp([all_data(i).name, ' did not converge']);
    end
end

% remove zeros
t(t==0) = [];
s_inf(t==0) = [];
s_one(t==0) = [];
s_freq(t==0) = [];

% keep only ones with 0 slack freq
% t(s_freq ~= 0) = []; 
% s_inf(s_freq ~= 0) = [];
% s_one(s_freq ~= 0) = [];
% s_freq(s_freq ~= 0) = [];

% sort in ascending
[ts, si] = sort(t, 'ascend');
s_inf_s = s_inf(si); 
s_freq_s = s_freq(si); 

% keep middle 5
Ng = numel(t); 
indm = floor(Ng/2) + 1;
tf = ts(indm - 2: indm + 2);
s_inf_f = s_inf_s(indm - 2: indm + 2);
s_freq_f = s_freq_s(indm - 2: indm + 2); 

fprintf('Average time : %f +/- % f (n = %d)\r', mean(tf)/60, std(tf)/60, numel(tf));
fprintf('Average max slack : %f +/- % f  (n = %d) \r', mean(s_inf_f), std(s_inf_f), numel(s_inf_f));
fprintf('Average slack freq : %f +/- % f  (n = %d) \r', mean(s_freq_f), std(s_freq_f), numel(s_freq_f));

save('SFTimingResults', 'tf', 's_inf_f')