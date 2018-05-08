function [traj, brkVal] = GenerateHueristicActTraj(gait, freq, dc, dl, NSAMP, TYPE)

NL = 4;

switch gait
    case 'TROT'
        % HACKED
        %         legphase = [0, 0, 0, 0, 0, 0, 0, 0]/(2*pi);     % fraction of period
        actphase = [0, 0, pi, pi, pi, pi, 0, 0]/(2*pi);     % fraction of period

    case 'PRONK'
         actphase = [0, 0, 0, 0, 0, 0, 0, 0]/(2*pi);     % fraction of period
%         actphase = [pi/2, 0, pi/2, pi, -pi/2, pi, -pi/2, 0]/(2*pi);
    case 'JUMP'
        actphase = [0, 0, pi/2, pi/2, 0, 0, pi/2, pi/2]/(2*pi);
    case 'BOUND' %check
        actphase = [pi/2, 0, pi/2, 0, -pi/2, pi, -pi/2, pi];
end

DIRECTIONS = [1, 1; -1, 1; -1, -1; 1, -1];
traj = zeros(NSAMP, numel(DIRECTIONS));
brkVal = cell(numel(DIRECTIONS),1);

for i = 1:NL
    
    lift_phi = actphase(2*i);
    swing_phi = actphase(2*i-1);
    
    switch TYPE
        case 1            % first type (varied impedance)
            
            [t, traj(:,2*i-1), brkVal{2*i-1}] = variable_dc_sine_fit(freq, dc, dl, DIRECTIONS(i,1), swing_phi, NSAMP, 'swing');
            [~, traj(:,2*i), brkVal{2*i}] = variable_dc_sine_fit(freq, dc, dl, DIRECTIONS(i,2), lift_phi, NSAMP, 'lift');
            
        case 2            % second type (varied impulse)
            
            [t, traj(:,2*i-1), brkVal{2*i-1}] = variable_dc_sine_fit2(freq, dc, dl, DIRECTIONS(i,1), swing_phi, NSAMP, 'swing');
            [~, traj(:,2*i), brkVal{2*i}] = variable_dc_sine_fit2(freq, dc, dl, DIRECTIONS(i,2), lift_phi, NSAMP, 'lift');
            
    end
    
end

traj = [t,traj];
end


function [tout, yy, brkVal] = variable_dc_sine_fit(f, dc, dl, direction, phase, NSAMP, state)
% frequency in Hz
% duty cycle is 0-100
% lift push is 0-100
% direction is orientation of signal (+/-)
% phase is in % cycle
% NSAMP is # of samples
% state determine if it's lift or swing

dc = dc/100;                % duty cycle
dl = dl/100;                % normal push
tn = linspace(0, 1, NSAMP)';         % one period of normalized time

tknots = [0, (1-dc)/4, 3*(1-dc)/4, 1-dc, 1-3*dc/4,  1-dc/4, 1];
qknots = direction*[-1, -1/2,       1/2,      1,   1/2,  -1/2, -1];
vknots = direction*[0,  2/(1-dc),  2/(1-dc), 0, -2/dc,  -2/dc, 0];
pps = hermite_spline(tknots, qknots, vknots);
swing = ppval(pps, tn);

tknotl = [0, 1/8, 1/4, 3/8, 1/2, 1];
qknotl = direction*[-dl, -dl+(1+dl)/2, 1, -dl+(1+dl)/2, -dl, -dl];
vknotl = direction*[0, 4*(1+dl), 0, -4*(1+dl), 0, 0];
ppl = hermite_spline(tknotl, qknotl, vknotl);
lift = ppval(ppl, tn);

%         figure(99); hold on;
%         plot(tn, swing);
%         plot(tn, lift);
%         plot(tknots, qknots, 'rx-');
%         plot(tknotl, qknotl, 'gx-');

swing = circshift(swing, phase*numel(tn));
lift = circshift(lift, phase*numel(tn));

% build time
T = 1/f;                            % period
tout = linspace(0, T, NSAMP)';      % time

switch state
    case 'lift'
        yy = lift;
        brkVal = [tknotl; qknotl; vknotl];
    case 'swing'
        yy = swing;
        brkVal = [tknots; qknots; vknots];
end

end

function [tout, yy, brkVal] = variable_dc_sine_fit2(f, dcs, dcl, direction, phase, NSAMP, state)
% frequency in Hz
% duty cycle swing is 0-100
% duty cycle lift is 0-100
% direction is orientation of signal (+/-)
% phase is in % cycle
% NSAMP is # of samples
% state determine if it's lift or swing

dcs = dcs/100;        % duty cycle swing
dcl = dcl/100;        % duty cycle lift
tn = linspace(0, 1, NSAMP)';         % one period of normalized time

tknots = [0, (1-dcs)/4, 3*(1-dcs)/4, 1-dcs, 1-3*dcs/4,  1-dcs/4, 1];
qknots = direction*[-1, -1/2,       1/2,      1,   1/2,  -1/2, -1];
vknots = direction*[0,  2/(1-dcs),  2/(1-dcs), 0, -2/dcs,  -2/dcs, 0];
pps = hermite_spline(tknots, qknots, vknots);
swing = ppval(pps, tn);

tknotl = [0, (1-dcl)/4, 3*(1-dcl)/4, 1-dcl, 1-3*dcl/4,  1-dcl/4, 1];
qknotl = direction*[-1, -1/2,       1/2,      1,   1/2,  -1/2, -1];
vknotl =direction*[0,  2/(1-dcl),  2/(1-dcl), 0, -2/dcl,  -2/dcl, 0];
ppl = hermite_spline(tknotl, qknotl, vknotl);
lift = ppval(ppl, tn);

swing = circshift(swing, phase*numel(tn));
lift = circshift(lift, (phase - 0.25)*numel(tn));

%         figure(99); hold on;
%         plot(tn, swing);
%         plot(tn, lift);
%         plot(tknots, qknots, 'rx-');
%         plot(tknotl, qknotl, 'gx-');


% build time
T = 1/f;                            % period
tout = linspace(0, T, NSAMP)';      % time

switch state
    case 'lift'
        yy = lift;
        brkVal = [tknotl; qknotl; vknotl];
    case 'swing'
        yy = swing;
        brkVal = [tknots; qknots; vknots];
end

end

function pp = hermite_spline(tknot, qknot, vknot)


N = numel(tknot) -1;
coeffs = zeros(4, N);
for i = 1:N
    t0 = tknot(i); t1 = tknot(i+1);
    q0 = qknot(i);  q1 = qknot(i+1);
    v0 = vknot(i);  v1 = vknot(i+1);
    
    h = t1 - t0;
    
    Ah = [2/h^3,  1/h^2 , -2/h^3,  1/h^2;
        -3/h^2, -2/h,  3/h^2, -1/h;
        0,    1,    0,  0;
        1,    0,    0,  0];
    
    coeffs(:,i) = Ah*[q0; v0; q1; v1];
    
end

pp = mkpp(tknot, coeffs');



end
