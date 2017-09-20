function A = randomSymmetricPositiveDefiniteMatrix(n) %#codegen
%RANDOMSYMMETRICPOSITIVEDEFINITEMATRIX Construct random SPD matrix.
%
% Syntax:
%   A = randomSymmetricPositiveDefiniteMatrix(n);
%
% Description:
%   Constructs a symmetric positive definite square matrix by making it
%   diagonally dominant and using A + A' to make it symmetric.
%
% Notes:
%   A + A' is computationally more efficient than A*A'.

% Copyright 2016-2017 Mikhail S. Jones

  % Set code generation settings
  coder.allowpcode('plain');
      
  % Initialize matrix
  A = rand(n);

  % Make matrix symmetric
  A = A + A.';

  % Make matrix positive definite
  A = A + n*eye(n);
end % randomSymmetricPositiveDefiniteMatrix
