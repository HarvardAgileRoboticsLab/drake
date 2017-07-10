clear; clc; close all; 
format shortg

sigma = 0.042e-2;  % g/mm^2

links = [10; 5];  % mm
output = [10; 10];  % mm

mlink = prod(links)*sigma        %g
moutput = prod(output)*sigma      %g

Ilink = (mlink/12)*(links(1)^2 + links(2)^2)
Ioutput = (moutput/12)*(output(1)^2 + output(2)^2)


%% flexures 

w= 5e-3; l = 200e-6; t = 12.5e-6; E = 2.5e9;

k = E*w*t^3/12/l
ks = k*1e3

b = 0.05*k/(2*pi*1);
bs = b*1e6;





