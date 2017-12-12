clear; clc; close all; 

trial_name = 'TrajOpt_MovingBody_SimpleSprings10';
save_dir = '~/Dropbox/CurrentWork/FrictionTrajOpt/MatFiles/TrajOptFiles/';

[hamr,xtraj,utraj,ctraj,btraj,...
    psitraj,etatraj,jltraj, kltraj, straj, ...
    z,F,info,infeasible_constraint_name] = SimpleHAMRVariationalTrajOpt(save_dir); 

save([save_dir, trial_name], 'xtraj', 'utraj', 'ctraj', 'btraj', 'psitraj', 'etatraj', ...
    'jltraj', 'kltraj', 'straj')