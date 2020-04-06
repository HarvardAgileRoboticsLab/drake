clear; clc; close all; 

urdf = fullfile(getDrakePath,'examples', 'HAMR-URDF', 'urdf', 'FLwithtwist.urdf');
twist_test = TimeSteppingRigidBodyManipulator(urdf, 1);
v = twist_test.constructVisualizer(); 
x0 = zeros(2*twist_test.getNumPositions, 1);
v.inspector(x0)


% xtraj = twist_test.simulate([0, 100], x0)
% tt = xtraj.getBreaks();
% xx = xtraj.eval(tt); 
% 
% plot(tt, xx(2,:))
% v.playback(xtraj)
% v.inspector(x0)