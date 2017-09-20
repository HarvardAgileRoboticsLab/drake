function A = randomPositiveDefiniteMatrix(n) %#codegen
%RANDOMPOSITIVEDEFINITEMATRIX Construct random positive definite matrix.
%
% Syntax:
%   A = randomPositiveDefiniteMatrix(n);
%
% Description:
%   Constructs a positive definite square matrix by making it diagonally
%   dominant.

% Copyright 2016-2017 Mikhail S. Jones

  % Set code generation settings
  coder.allowpcode('plain');
  
  % Initialize matrix
  A = rand(n);

  % Make matrix positive definite
  A = A + n*eye(n);
end % randomPositiveDefiniteMatrix
