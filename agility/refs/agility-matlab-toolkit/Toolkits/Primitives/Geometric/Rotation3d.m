%ROTATION3D Represents a rotation in three dimensional space.
%
% Syntax:
%   rotation = Rotation3d();
%   rotation = Rotation3d(eye(3));

% Copyright 2016-2017 Mikhail S. Jones

classdef Rotation3d %#codegen

  % PROTECTED PROPERTIES ==================================================
  properties (Access = protected)
    % Rotation matrix data
    data
  end % properties

  % PUBLIC METHODS ========================================================
  methods (Access = public)
    function obj = Rotation3d(data)
    %ROTATION3D Construct a rotation in three dimensions.

      % Set code generation settings
      coder.allowpcode('plain');
      
      % Check for no input arguments
      if nargin == 0
        % Set zero rotation data
        obj.data = eye(3);
      else
        % Validate rotation data attributes
        validateattributes(data, ...
          {'double', 'sym'}, ...
          {'real', 'size', [3, 3]}, ...
          'Rotation3d', ...
          'data');
        
        % Set rotation data
        obj.data = data;
      end % if
    end % Rotation3d
    
    function disp(obj)
    %DISP Display object value.
      disp(obj.getValue);
    end % disp
    
    % MATH OPERATOR METHODS -----------------------------------------------
    function a = uplus(a)
    %UPLUS Unary addition.
    end % uplus

    function a = uminus(a)
    %UMINUS Unary subtraction.
      a.data = -a.data;
    end % uminus
    
    function c = mtimes(a, b)
    %MTIMES Matrix multiply.

      % Check datatypes
      if isa(a, 'Rotation3d') && isa(b, 'Rotation3d')
        % Rotation = Rotation * Rotation
        c = Rotation3d(a.getValue*b.getValue);
      elseif isa(a, 'Rotation3d') && isa(b, 'Vector3d')
        % Vector = Rotation * Vector
        c = Vector3d(a.getValue*b.getValue);
      else
        % Throw error
        error('Unsupported multiplication.');
      end % if
    end % mtimes

    function c = mrdivide(a, b)
    %MRDIVIDE Matrix right handed division.
    %
    % Note:
    %   - Uses inv(R) == R' assumption

      % Check datatypes
      if isa(a, 'Rotation3d') && isa(b, 'Rotation3d')
        % Rotation = Rotation / Rotation
        c = Rotation3d(a.getValue*b.getValue.');
      else
        % Throw error
        error('Unsupported division.');
      end % if
    end % mrdivide

    function c = mldivide(a, b)
    %MLDIVIDE Matrix left handed division.
    %
    % Note:
    %   - Uses inv(R) == R' assumption
    
      % Check datatypes
      if isa(a, 'Rotation3d') && isa(b, 'Rotation3d')
        % Rotation = Rotation \ Rotation
        c = Rotation3d(a.getValue.'*b.getValue);
      elseif isa(a, 'Rotation3d') && isa(b, 'Vector3d')
        % Vector = Rotation \ Vector
        c = Vector3d(a.getValue.'*b.getValue);
      else
        % Throw error
        error('Unsupported division.');
      end % if
    end % mldivide
    
    function a = inverse(a)
    %INVERSE Matrix inverse.
    %
    % Note:
    %   - Uses inv(R) == R' assumption
      a.data = a.data.';
    end % inverse
    
    function a = transpose(a)
    %TRANSPOSE Matrix transpose.
      a.data = a.data.';
    end % transpose

    function t = trace(obj)
    %TRACE Sum of the diagonal components.
    
      % Extract rotation matrix
      R = obj.data;
      
      % Compute trace
      t = R(1,1) + R(2,2) + R(3,3);
    end % trace
    
    function v = skew(obj)
    %SKEW Skew operation.

      % Extract rotation matrix
      R = obj.data;

      % Compute skew vector
      v = [R(3,2) - R(2,3); R(1,3) - R(3,1); R(2,1) - R(1,2)]/2;
    end % skew

    % ROTATION OPERATOR METHODS -------------------------------------------
    function obj = rotX(obj, rx)
    %ROTX Rotate around X-axis.
    %
    % Syntax:
    %   Rotation().rotX(rx);

      % Validate rotation angle attributes
      validateattributes(rx, ...
        {'double', 'sym'}, ...
        {'real', 'scalar'}, ...
        'Rotation3d.rotX', ...
        'rx');

      % Compute intermediate variables
      c = cos(rx);
      s = sin(rx);

      % Construct rotation matrix
      Rx = [1, 0, 0; 0, c, -s; 0, s, c];

      % Update rotation object
      obj.data = obj.data*Rx;
    end % rotX

    function obj = rotY(obj, ry)
    %ROTY Rotate around Y-axis.
    %
    % Syntax:
    %   Rotation().rotY(ry);

      % Validate rotation angle attributes
      validateattributes(ry, ...
        {'double', 'sym'}, ...
        {'real', 'scalar'}, ...
        'Rotation3d.rotY', ...
        'ry');

      % Compute intermediate variables
      c = cos(ry);
      s = sin(ry);

      % Construct rotation matrix
      Ry = [c, 0, s; 0, 1, 0; -s, 0, c];

      % Update rotation object
      obj.data = obj.data*Ry;
    end % rotY

    function obj = rotZ(obj, rz)
    %ROTZ Rotate around Z-axis.
    %
    % Syntax:
    %   Rotation().rotZ(rz);

      % Validate rotation angle attributes
      validateattributes(rz, ...
        {'double', 'sym'}, ...
        {'real', 'scalar'}, ...
        'Rotation3d.rotZ', ...
        'rz');

      % Compute intermediate variables
      c = cos(rz);
      s = sin(rz);

      % Construct rotation matrix
      Rz = [c, -s, 0; s, c, 0; 0, 0, 1];

      % Update rotation object
      obj.data = obj.data*Rz;
    end % rotZ

    function obj = rotZYX(obj, r)
    %ROTZYX Rotate using ZYX Euler angles.
    %
    % Syntax:
    %   Rotation().rotZYX(r);

      % Validate rotation angle attributes
      validateattributes(r, ...
        {'double', 'sym'}, ...
        {'real', 'numel', 3}, ...
        'Rotation3d.rotZYX', ...
        'r');

      % Update rotation object
      obj = obj.rotZ(r(1)).rotY(r(2)).rotX(r(3));
    end % rotZYX

    function obj = rotZYZ(obj, r)
    %ROTZYZ Rotate using ZYZ Euler angles.
    %
    % Syntax:
    %   Rotation().rotZYZ(r);

      % Validate rotation angle attributes
      validateattributes(r, ...
        {'double', 'sym'}, ...
        {'real', 'numel', 3}, ...
        'Rotation3d.rotZYZ', ...
        'r');

      % Update rotation object
      obj = obj.rotZ(r(1)).rotY(r(2)).rotZ(r(3));
    end % rotZYZ
    
    function obj = rotAxis(obj, axis, angle)
    %ROTAXIS Rotate about arbitrary axis.
    %
    % Syntax:
    %   Rotation().rotAxis(axis, angle);

      % Validate axis vector attributes
      validateattributes(axis, ...
        {'Vector'}, ...
        {'scalar'}, ...
        'Rotation3d.rotAxis', ...
        'axis');
      
      % Validate rotation angle attributes
      validateattributes(angle, ...
        {'double', 'sym'}, ...
        {'real', 'scalar'}, ...
        'Rotation3d.rotAxis', ...
        'angle');

      % Compute intermediate variables
      c = cos(angle);
      s = sin(angle);

      % Normalize axis
      axis = axis.normalize;
      
      % Break out vector components
      x = axis.getX;
      y = axis.getY;
      z = axis.getZ;      
      
      % Construct rotation matrix
      R = [c + x^2*(1 - c), x*y*(1 - c) - z*s, x*z*(1 - c) + y*s;
        y*x*(1 - c) + z*s, c + y^2*(1 - c), y*z*(1 - c) - x*s;
        z*x*(1 - c) - y*s, z*y*(1 - c) + x*s, c + z^2*(1 - c)];

      % Update rotation object
      obj.data = obj.data*R;
    end % rotAxis
    
    function obj = removePitchAndRoll(obj)
    %REMOVEPITCHANDROLL Remove pitch and roll rotations from matrix.
      
      % Extract the rotation matrix
      R = obj.data;
      
      % Extract rotation matrix components responsible for yaw
      R11 = R(1,1);
      R21 = R(2,1);
      
      % Compute the norm
      n = sqrt(R11^2 + R21^2);
      
      % Check if norm is zero to prevent dividing by zero
      if n == 0
        % Identity rotation matrix
        obj.data = eye(3);
      else
        % Compute scaled components of rotation matrix
        R11 = R11/n;
        R21 = R21/n;
        
        % Construct projected rotation matrix
        obj.data = [...
          R11, -R21, 0; ...
          R21, R11, 0; ...
          0, 0, 1];
      end % if
    end % removePitchAndRoll
  end % methods
  
  % GET/SET METHODS =======================================================
  methods
    % Rotation ------------------------------------------------------------
    function value = getValue(obj)
    %GETVALUE Get value.
      value = obj.data;
    end % getValue
    
    % Quaternion ----------------------------------------------------------
    function q = getQuaternion(obj)
    %GETQUATERNION Convert rotation matrix to quaternion.
    %
    % Notes:
    %   - https://en.wikipedia.org/wiki/Rotation_matrix#Quaternion
    %   - http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
    
      % Extract the rotation matrix
      R = obj.data;

      % Compute the trace of the rotation matrix
      trace = obj.trace;
  
      % Check which major diagonal element has the greatest value if the
      % trace is less than or equal to zero
      if trace > 0
        % Compute the quaternion
        S = 2*sqrt(trace + 1);
        qw = 0.25*S;
        qx = (R(3,2) - R(2,3))/S;
        qy = (R(1,3) - R(3,1))/S; 
        qz = (R(2,1) - R(1,2))/S;
        
      elseif R(1,1) > R(2,2) && R(1,1) > R(3,3)
        % Compute the quaternion
        S = 2*sqrt(1 + R(1,1) - R(2,2) - R(3,3));
        qw = (R(3,2) - R(2,3))/S;
        qx = 0.25*S;
        qy = (R(1,2) + R(2,1))/S; 
        qz = (R(1,3) + R(3,1))/S;
        
      elseif R(2,2) > R(3,3)
        % Compute the quaternion
        S = 2*sqrt(1 + R(2,2) - R(1,1) - R(3,3));
        qw = (R(1,3) - R(3,1))/S;
        qx = (R(1,2) + R(2,1))/S; 
        qy = 0.25*S;
        qz = (R(2,3) + R(3,2))/S;
        
      else
        % Compute the quaternion
        S = 2*sqrt(1.0 + R(3,3) - R(1,1) - R(2,2));
        qw = (R(2,1) - R(1,2))/S;
        qx = (R(1,3) + R(3,1))/S;
        qy = (R(2,3) + R(3,2))/S;
        qz = 0.25*S;
      end % if
      
      % Construct quaternion
      q = Quaternion([qw; qx; qy; qz]);
    end % getQuaternion
    
    % Euler Angles --------------------------------------------------------
    function eulerAngles = getEulerAngles(obj)
    %GETEULERANGLES Convert the rotation matrix to Euler angles.
      
      % Construct vector of Euler angles
      eulerAngles = [...
        obj.getYaw;
      	obj.getPitch;
        obj.getRoll];
    end % getEulerAngles
    
    function yaw = getYaw(obj)
    %GETYAW Compute the yaw.
    
      % Extract the rotation matrix
      R = obj.data;
      
      % Compute yaw
      yaw = atan2(R(2,1), R(1,1));
    end % getYaw
    
    function pitch = getPitch(obj)
    %GETPITCH Compute the pitch.
    
      % Extract the rotation matrix
      R = obj.data;
      
      % Compute pitch
      pitch = -asin(R(3,1));
    end % getPitch
    
    function roll = getRoll(obj)
    %GETROLL Compute the roll.
    
      % Extract the rotation matrix
      R = obj.data;
      
      % Compute roll
      roll = atan2(R(3,2), R(3,3));
    end % getRoll
    
    % Axis-Angle ----------------------------------------------------------
    function [axis, angle] = getAxisAngle(obj)
    %GETAXISANGLE Convert the rotation matrix to axis-angle form.
    %
    % Notes:
    %   http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToAngle/
    
      % Extract the rotation matrix
      R = obj.data;
      
      % Compute the trace of the rotation matrix
      trace = obj.trace;
      
      % Compute the angle
      angle = acos((trace - 1)/2);
      
      % Compute the axis
      axis = [...
        R(3,2) - R(2,3); ...
        R(1,3) - R(3,1); ...
        R(2,1) - R(1,2)]/(2*sin(angle)); 
    end % getAxisAngle
  end % methods

  % STATIC METHODS ========================================================
  methods (Static = true)
    function obj = rand
    %RAND Construct uniformly distributed pseudorandom rotation.
      obj = rotZYX(Rotation3d, pi*rand(3,1));
    end % rand
    
    function obj = randn
    %RANDN Construct normally distributed pseudorandom rotation.
      obj = rotZYX(Rotation3d, pi*randn(3,1));
    end % randn
  end % methods
end % classdef
