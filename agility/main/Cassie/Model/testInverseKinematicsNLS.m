function [q, status] = testInverseKinematicsNLS(leftFootTransform, rightFootTransform)
%TESTINVERSEKINEMATICSNLS Test inverse kinematics NLS solver.
%
% Syntax:
%   [q, status] = testInverseKinematicsNLS
%
% Description:
%   Formulates the inverse kinematics problem as a nonlinear least squares
%   optimization problem with joint limit constraints.

% Copyright 2016-2017 Mikhail S. Jones

  % Get the motor position limits
  lowerBound = CassieParameters.motorPositionLowerLimit;
  upperBound = CassieParameters.motorPositionUpperLimit;
  
  % Construct the inverse kinematics function
  inverseKinematicsFunction = InverseKinematicsFunction;
  
  % Construct the nonlinear least squares problem
  problem = NonlinearLeastSquaresProblem(...
    inverseKinematicsFunction, lowerBound, upperBound);

  % Construct the nonlinear least squares solver
  solver = NonlinearLeastSquaresSolver(problem);
  solver.setIterationLimit(100);
  solver.setGTolerance(eps);
  solver.setXTolerance(eps);

  % Check the number of input arguments
  if nargin == 0
    % Define the world transform
    worldTransform = Transform3d;
  
    % Define the desired toe transforms
    leftFootTransform = worldTransform.translate([0; 0.135; -1]);
    rightFootTransform = worldTransform.translate([0; -0.135; -1]);
  end % if
    
  % Update the desired foot transforms in the inverse kinematics problem
  inverseKinematicsFunction.setDesiredFootTransforms(...
    leftFootTransform, rightFootTransform);
    
  % Solve the problem
  solver.solve;
  
  % Get the solver status
  status = solver.getStatus;
  
  % Get the solver solution
  q = solver.getSolution;
end % testInverseKinematicsNLS