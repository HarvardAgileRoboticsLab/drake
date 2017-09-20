function [q, dq, ddq] = cubicInterpolation(t, ts, qs, dqs)
%CUBICINTERPOLATION Cubic interpolation between vectors.

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
        1, t0, t0^2, t0^3;
        0, 1, 2*t0, 3*t0^2;
        1, tf, tf^2, tf^3;
        0, 1, 2*tf, 3*tf^2];

      % Construct system of equations vector
      b = [...
        qs(i);
        dqs(i);
        qs(i+1);
        dqs(i+1)];

      % Solve system of equations for coefficients
      c = A\b;

      % Compute interpolated position
      q = c(1) + c(2)*t + c(3)*t^2 + c(4)*t^3;

      % Compute interpolated velocity
      dq = c(2) + 2*c(3)*t + 3*c(4)*t^2;

      % Compute interpolated acceleration
      ddq = 2*c(3) + 6*c(4)*t;
    end % if
  end % for
  
  % Set acceleration to zero at end points
  if t <= ts(1) || t >= ts(end)
    ddq = 0;
  end % if
end % cubicInterpolation
