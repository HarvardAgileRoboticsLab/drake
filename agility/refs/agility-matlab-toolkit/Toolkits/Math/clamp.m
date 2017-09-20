function b = clamp(a, limit_1, limit_2) %#codegen
%CLAMP Clamp value between two bounds.
%
% Syntax:
%   b = clamp(a, limit_1, limit_2);
%
% Todo:
%   - Rename to clip

% Copyright 2015-2017 Mikhail S. Jones

  % Set code generation settings
  coder.allowpcode('plain');
      
  % Validate input argument a attributes
  validateattributes(a, ...
    {'double'}, ...
    {'real', 'vector'}, ...
    'clamp', ...
    'a');
  
  % Validate input argument limit_1 attributes
  validateattributes(limit_1, ...
    {'double'}, ...
    {'real', 'scalar'}, ...
    'clamp', ...
    'limit_1');
  
  % Validate input argument limit_2 attributes
  validateattributes(limit_2, ...
    {'double'}, ...
    {'real', 'scalar'}, ...
    'clamp', ...
    'limit_2');
  
  % Determine lower and upper limits
  lowerLimit = min(limit_1, limit_2);
  upperLimit = max(limit_1, limit_2);

  % Clamp value between limits
  b = max(min(a, upperLimit), lowerLimit);
end % clamp
