%QUATERNION Represents a quaternion.
%
% Syntax:
%   quaternion = Quaternion();
%   quaternion = Quaternion([w; x; y; z]);
%
% TODO:
%   - Add quaternion differentiation method diff(qVec, dt)
%   - Add mean method
%   - Add quaternion log and exp methods

% Copyright 2016-2017 Mikhail S. Jones

classdef Quaternion %#codegen

  % PROTECTED PROPERTIES ==================================================
  properties (Access = protected)
    % Quaternion data
    data
  end % properties

  % PUBLIC METHODS ========================================================
  methods (Access = public)
    function obj = Quaternion(data)
    %QUATERNION Construct a quaternion.

      % Set code generation settings
      coder.allowpcode('plain');
      
      % Check for no input arguments
      if nargin == 0
        % Set zero quaternion data
        obj.data = [1; 0; 0; 0];
      else
        % Validate quaternion data attributes
        validateattributes(data, ...
          {'double', 'sym'}, ...
          {'real', 'size', [4, 1]}, ...
          'Quaternion', ...
          'data');
        
        % Set quaternion data
        obj.data = data;
      end % if
    end % Quaternion
    
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

    function a = plus(a, b)
    %PLUS Addition.

      % Check input argument datatypes
      if isa(a, 'Quaternion') && isa(b, 'Quaternion')
        % Element wise vector addition
        a.data = a.data + b.data;
      else
        error('Unsupported addition.');
      end % if
    end % plus

    function a = minus(a, b)
    %MINUS Subtraction.

      % Check input argument datatypes
      if isa(a, 'Quaternion') && isa(b, 'Quaternion')
        % Element wise vector subtraction
        a.data = a.data - b.data;
      else
        error('Unsupported subtraction.');
      end % if
    end % minus
    
    function c = mtimes(a, b)
    %MTIMES Matrix multiply.

      % Check input argument datatypes
      if isa(a, 'Quaternion') && isa(b, 'Quaternion')
        % Get local variables
        aw = a.getW; ax = a.getX; ay = a.getY; az = a.getZ;
        bw = b.getW; bx = b.getX; by = b.getY; bz = b.getZ;

        % Quaternion = Quaternion * Quaternion
        cw = bw*aw - bx*ax - by*ay - bz*az;
        cx = bw*ax + bx*aw - by*az + bz*ay;
        cy = bw*ay + bx*az + by*aw - bz*ax;
        cz = bw*az - bx*ay + by*ax + bz*aw;

        % Construct quaternion
        c = Quaternion([cw; cx; cy; cz]);
        
      elseif isa(a, 'double') && numel(a) == 1
        % Scalar multiplication
        c = Quaternion(a*b.data);        
      elseif isa(b, 'double') && numel(b) == 1
        % Scalar multiplication
        c = Quaternion(a.data*b);        
      else
        % Throw error
        error('Unsupported multiplication.');
      end % if
    end % mtimes
    
    function a = mrdivide(a, b)
    %MRDIVIDE Matrix right divide.

      % Check input argument datatypes
      if isa(b, 'double') && numel(b) == 1
        % Scalar division
        a.data = a.data/b;
      else
        error('Unsupported division.');
      end % if
    end % mrdivide

    function n = norm(obj)
    %NORM Compute norm of quaternion.
      n = sum(obj.data.^2);
    end % norm
    
    function obj = normalize(obj)
    %NORMALIZE Normalize quaternion.
      obj.data = obj.data/sqrt(obj.norm);
    end % normalize
    
    function a = conj(a)
    %CONJ Compute conjugate of quaternion.
      a.data(2:4) = -a.data(2:4);
    end % conj
  end % methods

  % GET/SET METHODS =======================================================
  methods
    % Quaternion ----------------------------------------------------------
    function value = getValue(obj)
    %GETVALUE Get value.
      value = obj.data;
    end % getValue
    
    function w = getW(obj)
    %GETW Get method for w-coordinate.
      w = obj.data(1);
    end % getW
    
    function x = getX(obj)
    %GETX Get method for x-coordinate.
      x = obj.data(2);
    end % getX

    function y = getY(obj)
    %GETY Get method for y coordinate.
      y = obj.data(3);
    end % getY

    function z = getZ(obj)
    %GETZ Get method for z-coordinate.
      z = obj.data(4);
    end % getZ
    
    % Rotation ------------------------------------------------------------
    function rotation = getRotation(obj)
    %GETROTATION Convert quaternion to rotation matrix.
      
      % Get local variables
      qw = obj.getW;
      qx = obj.getX;
      qy = obj.getY;
      qz = obj.getZ;
      
      % Compose homogeneous rotation matrix
      rotationMatrix = [...
        qw^2 + qx^2 - qy^2 - qz^2, 2*(qx*qy - qw*qz), 2*(qw*qy + qx*qz);
        2*(qx*qy + qw*qz), qw^2 - qx^2 + qy^2 - qz^2, 2*(qy*qz - qw*qx);
        2*(qx*qz - qw*qy), 2*(qw*qx + qy*qz), qw^2 - qx^2 - qy^2 + qz^2];
      
      % Construct rotation object
      rotation = Rotation3d(rotationMatrix);
    end % getRotation
    
    % Euler Angles --------------------------------------------------------
    function eulerAngles = getEulerAngles(obj)
    %GETEULERANGLES Convert quaternion to Euler angles.
      
      % Construct vector of Euler angles
      eulerAngles = [...
        obj.getYaw;
      	obj.getPitch;
        obj.getRoll];
    end % getEulerAngles
    
    function yaw = getYaw(obj)
    %GETYAW Compute yaw.
    %
    % Notes:
    %   - Singularity when pitch equals plus or minus pi/2
      
      % Get local variables
      qw = obj.getW;
      qx = obj.getX;
      qy = obj.getY;
      qz = obj.getZ;
      
      % Compute yaw
      yaw = atan2(2*(qx*qy + qw*qz),  qw^2 + qx^2 - qy^2 - qz^2);
    end % getYaw
    
    function pitch = getPitch(obj)
    %GETPITCH Compute pitch.
    
      % Get local variables
      qw = obj.getW;
      qx = obj.getX;
      qy = obj.getY;
      qz = obj.getZ;
      
      % Compute pitch
      pitch = asin(-2*(qx*qz - qw*qy));
    end % getPitch
    
    function roll = getRoll(obj)
    %GETROLL Compute roll.
    %
    % Notes:
    %   - Singularity when pitch equals plus or minus pi/2
    
      % Get local variables
      qw = obj.getW;
      qx = obj.getX;
      qy = obj.getY;
      qz = obj.getZ;
      
      % Compute roll
      roll = atan2(2*(qw*qx + qy*qz), qw^2 - qx^2 - qy^2 + qz^2);
    end % getRoll
  end % methods

  % STATIC METHODS ========================================================
  methods (Static = true)
    function obj = rand
    %RAND Construct random quaternion object.
      obj = Quaternion(rand(4,1)).normalize;
    end % rand
  end % methods
end % classdef
