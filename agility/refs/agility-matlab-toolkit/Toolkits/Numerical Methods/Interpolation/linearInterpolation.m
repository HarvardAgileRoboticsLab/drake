function [q, dq, ddq] = linearInterpolation(t, ts, qs)
%LINEARINTERPOLATION Linear interpolation between vectors.

% Copyright 2016-2017 Mikhail S. Jones
  
  % Validate time vector
  validateattributes(ts, ...
    {'double'}, ...
    {'vector', 'increasing'}, ...
    'linearInterpolation', ...
    'ts');
  
  % Validate position vector
  validateattributes(qs, ...
    {'double'}, ...
    {'vector', 'numel', numel(ts)}, ...
    'linearInterpolation', ...
    'qs');
  
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
        1, t0;
        1, tf];

      % Construct system of equations vector
      b = [...
        qs(i);
        qs(i+1)];

      % Solve system of equations for coefficients
      c = A\b;

      % Compute interpolated position
      q = c(1) + c(2)*t;

      % Compute interpolated velocity
      dq = c(2);

      % Compute interpolated acceleration
      ddq = 0;
    end % if
  end % for
end % linearInterpolation
