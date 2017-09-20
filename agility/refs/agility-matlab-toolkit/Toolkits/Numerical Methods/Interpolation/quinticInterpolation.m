function [q, dq, ddq] = quinticInterpolation(t, ts, qs, dqs, ddqs)  
%QUINTICINTERPOLATION Quintic interpolation between vectors.

% Copyright 2016-2017 Mikhail S. Jones
  
  % Validate time vector
  validateattributes(ts, ...
    {'double'}, ...
    {'vector', 'increasing'}, ...
    'cubicInterpolation', ...
    'ts');
  
  % Validate position vector
  validateattributes(qs, ...
    {'double'}, ...
    {'vector', 'numel', numel(ts)}, ...
    'cubicInterpolation', ...
    'qs');
  
  % Check number of input arguments
  if nargin < 4
    % Set default velocity vector if not supplied
    dqs = zeros(size(qs));
  end % if
  
  % Validate velocity vector
  validateattributes(dqs, ...
    {'double'}, ...
    {'vector', 'numel', numel(ts)}, ...
    'cubicInterpolation', ...
    'dqs');
  
  % Check number of input arguments
  if nargin < 5
    % Set default acceleration vector if not supplied
    ddqs = zeros(size(qs));
  end % if

  % Validate acceleration vector
  validateattributes(ddqs, ...
    {'double'}, ...
    {'vector', 'numel', numel(ts)}, ...
    'cubicInterpolation', ...
    'ddqs');
  
  % Limit interpolation point to be within trajectory domain
  t = median([ts(1), t, ts(end)]);

  % Initialize outputs
  q = qs(1); dq = 0; ddq = 0;

  % Loop through intervals
  for i = 1:numel(ts)-1
    % Find interval containing current point
    if t >= ts(i) && t <= ts(i+1)
      % Get initial and final interval times
      t0 = ts(i);
      tf = ts(i+1);

      % Construct system of equations matrix
      A = [...
        1, t0, t0^2, t0^3, t0^4, t0^5;
        0, 1, 2*t0, 3*t0^2, 4*t0^3, 5*t0^4;
        0, 0, 2, 6*t0, 12*t0^2, 20*t0^3;
        1, tf, tf^2, tf^3, tf^4, tf^5;
        0, 1, 2*tf, 3*tf^2, 4*tf^3, 5*tf^4;
        0, 0, 2, 6*tf, 12*tf^2, 20*tf^3];

      % Construct system of equations vector
      b = [...
        qs(i);
        dqs(i);
        ddqs(i);
        qs(i+1);
        dqs(i+1);
        ddqs(i+1)];

      % Solve system of equations for coefficients
      c = A\b;

      % Compute interpolated position
      q = c(1) + c(2)*t + c(3)*t^2 + c(4)*t^3 + c(5)*t^4 + c(6)*t^5;

      % Compute interpolated velocity
      dq = c(2) + 2*c(3)*t + 3*c(4)*t^2 + 4*c(5)*t^3 + 5*c(6)*t^4;

      % Compute interpolated acceleration
      ddq = 2*c(3) + 6*c(4)*t + 12*c(5)*t^2 + 20*c(6)*t^3;
    end % if
  end % for
end % quinticInterpolation