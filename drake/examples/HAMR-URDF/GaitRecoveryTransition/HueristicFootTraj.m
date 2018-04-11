function traj = HueristicFootTraj(gait, freq, dc, NPTS)

NL = 4;

switch gait
    case 'TROT'
        legphase = [-pi/2, 0, ...
            pi/2, pi, ...
            pi/2, pi, ... 
            -pi/2, 0];
    case 'PRONK' %check
        legphase = [pi/2, 0, pi/2, pi, -pi/2, pi, -pi/2, 0];
    case 'PACE' %check
        legphase = [pi/2, 0, pi/2, pi, pi/2, 0, pi/2, pi];
    case 'BOUND' %check
        legphase = [pi/2, 0, pi/2, 0, -pi/2, pi, -pi/2, pi];
end

swing_dc = dc;
lift_dc = 50;
swing_off = 0;
lift_off = 0; %legRad*cosd(legang)*(cosd(legamp)-1);
%legang = 30; %deg
%legamp = 10; %deg
%legRad = 10;

for i = 1:NL
    
    lift_phi = legphase(2*i-1);
    swing_phi = legphase(2*i);
    
    [t, traj(:,2*i-1)] = variable_dc_sine(freq, swing_dc, swing_phi, swing_off, NPTS);
    [~, traj(:,2*i)] = variable_dc_sine(freq, lift_dc, lift_phi, lift_off, NPTS);
    
end

%Normalizing to 1
traj = bsxfun(@rdivide, traj, max(traj));



traj = [t,traj];
end


function [t, y] = variable_dc_sine(f, dc, phase, offset, NPTS)
% frequency in Hz
% duty cycle is 0-100
% phase in radians
% offset (amplitude is 1);

N = 5;              % number of harmonics
dc = dc/100;        % duty cycle
T = 1/f;            % period
t = linspace(0, T, NPTS)';    % one period of time
t_cutoff1 = dc*T;   % cutoff time
sind= find(t >= t_cutoff1, 1, 'first'); %cutoff index

f_high = 0.5/(dc*T);        % frequency of signal above 0
f_low = 0.5/((1-dc)*T);     % frequency of signal below 0

T_high = 0.5/f_high;        % period of high signal
T_low = 0.5/f_low;          % period of low signal

% determine phase shift between high and low signals
dt = rem(T_high, T_low);
phi = -pi*dt/T_low;

%determine sign of low signal
ndiv = floor(T_high/T_low);
if mod(ndiv, 2) == 0
    sgn = -1;
else
    sgn = 1;
end

% Construct signal
x_high = cos(2*pi*f_high*t);
x_low = sgn*cos(2*pi*f_low*t + phi);

% truncate signal and time based on phase
dt = (T*(1 - 2*dc))/4;
if dt <= 0
    dt = T + dt;
end

sine = [x_high(1:sind); x_low(sind+1:end)];
t = [t; T+t(2:end)];
sine = [sine; sine(2:end)];

indl = find(t>=dt, 1, 'first');
indh = find(t>=dt+T, 1, 'first');
t = t(indl:indh)-t(indl);
sine = sine(indl:indh);

%Phase shift to make ensure
%phase=pi is mid-swing
%phase=0, 2*pi is mid-stance
zsine = crossing(sine);
zzsine = zsine((sine(zsine)-sine(zsine+1))<=0);
sine = circshift(sine, round(1.5*length(sine)-zzsine)-1);

%Additional phase, offset
sine = circshift(sine, round((phase*numel(t))/(2*pi))) + offset;


% fourier fit
% fourierfit = fit(t, sine, 'fourier5', opt);
% build regressors
f_fit = f*(1:N);
reg_fit = [0*t + 1, sin(2*pi*t*f_fit), cos(2*pi*t*f_fit)];
coeffs = reg_fit\sine;

% extract coefficients
% C0 = coeffs(1);             % constant
% As = coeffs(2:2:11);        % sine gains
% Bs = coeffs(3:2:12);        % cos gains
y = reg_fit*coeffs;

% figure(1); hold on;
% plot(repmat(sine, 1, 1));
% plot(repmat(y, 1, 1));


end
