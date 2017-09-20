function M = skewSymmetricMatrix(v) %#codegen
%SKEWSYMMETRICMATRIX Construct skew symmetric matrix.
%
% Syntax:
%   M = skewSymmetricMatrix(v);
%
% Description:
%   Constructs a skew symmetric square matrix.

% Copyright 2016-2017 Mikhail S. Jones

  % Set code generation settings
  coder.allowpcode('plain');
      
  % Validate vector attributes
  validateattributes(v, ...
    {'double', 'sym'}, ...
    {'numel', 3}, ...
    'skewSymmetricMatrix', ...
    'v');
  
  % Make skew symmetric matrix
  M = [...
    0, -v(3), v(2);
    v(3), 0, -v(1);
    -v(2), v(1), 0];
end % skewSymmetricMatrix
