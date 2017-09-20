%VECTOR3D Represents a vector in three dimensional space.
%
% Syntax:
%   v = Vector3d();
%   v = Vector3d([x; y; z]);

% Copyright 2016-2017 Mikhail S. Jones

classdef Vector3d %#codegen

  % PROTECTED PROPERTIES ==================================================
  properties (Access = protected)
    % Vector data
    data
  end % properties

  % PUBLIC METHODS ========================================================
  methods (Access = public)
    % CONSTRUCTOR ---------------------------------------------------------
    function obj = Vector3d(data)
    %VECTOR3D Construct a vector in three dimensional space.

      % Set code generation settings
      coder.allowpcode('plain');
      
      % Check for no input argument constructor
      if nargin == 0
        % Set zero vector data
        obj.data = zeros(3,1);
        
      else
        % Validate vector data
        validateattributes(data, ...
          {'double', 'sym'}, ...
          {'real', 'size', [3, 1]}, ...
          'Vector3d', ...
          'data');
        
        % Set vector data
        obj.data = data;
      end % if
    end % Vector3d
    
    function plot(obj, origin)
    %PLOT Plot vector.
    
      % Check for no input arguments
      if nargin <= 1
        origin = Transform3d;
      end % if
      
      % Scale of frame representation
      scale = 0.001;
      
      % Plot vector
      v = quiver3(0, 0, 0, scale*obj.getX, scale*obj.getY, scale*obj.getZ, ...
        'Color', 'm', ...
        'LineWidth', 3);
      
      % Construct handle graphics transform object
      hg = hgtransform;

      % Set plots to be children of handle graphics transform
      set(v, 'Parent', hg);
      
      % Set handle graphics object transformation matrix
      set(hg, 'Matrix', origin.getValue);
    end % plot
    
    function disp(obj)
    %DISP Display object value.
      disp(obj.getValue);
    end % disp

    % MATH OPERATORS ------------------------------------------------------
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
      if isa(a, 'Vector3d') && isa(b, 'Vector3d')
        % Element wise vector addition
        a.data = a.data + b.data;
      else
        error('Unsupported addition.');
      end % if
    end % plus

    function a = minus(a, b)
    %MINUS Subtraction.

      % Check input argument datatypes
      if isa(a, 'Vector3d') && isa(b, 'Vector3d')
        % Element wise vector subtraction
        a.data = a.data - b.data;
      else
        error('Unsupported subtraction.');
      end % if
    end % minus

    function c = mtimes(a, b)
    %MTIMES Matrix multiply.

      % Check input arguument datatypes
      if isa(a, 'double') && numel(a) == 1
        % Scalar multiplication
        c = Vector3d(a*b.data);
      elseif isa(b, 'double') && numel(b) == 1
        % Scalar multiplication
        c = Vector3d(a.data*b);
      else
        error('Unsupported multiplication.');
      end % if
    end % mtimes
    
    function c = times(a, b)
    %TIMES Array multiply.

      % Check input argument datatypes
      if isa(a, 'double') && numel(a) == 1
        % Scalar multiplication
        c = Vector3d(a*b.data);
      elseif isa(b, 'double') && numel(b) == 1
        % Sclar multiplication
        c = Vector3d(a.data*b);
      elseif isa(a, 'Vector3d') && isa(b, 'Vector3d')
        % Element wise vector multiplication
        c = Vector3d(a.data.*b.data);
      else
        error('Unsupported multiplication.');
      end % if
    end % times

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
    
    function a = rdivide(a, b)
    %RDIVIDE Array right divide.

      % Check input argument datatypes
      if isa(b, 'double') && numel(b) == 1
        % Scalar division
        a.data = a.data/b;
      elseif isa(a, 'Vector3d') && isa(b, 'Vector3d')
        % Element wise vector division
        a.data = a.data./b.data;
      else
        error('Unsupported division.');
      end % if
    end % rdivide

    function s = dot(a, b)
    %DOT Vector dot product.
    
      if isa(a, 'Vector3d') && isa(b, 'Vector3d')
        % Retun the scalar product of the vectors
        s = a.data.'*b.data;
      else
        error('Unsupported dot product.');
      end % if
    end % dot
    
    function n = norm(obj)
    %NORM Norm of vector.
      n = sqrt(obj.data.'*obj.data);
    end % norm
    
    function obj = normalize(obj)
    %NORMALIZE Normalize vector.
      obj.data = obj.data/norm(obj);
    end % normalize

    function S = skew(obj)
    %SKEW Skew operation.
      S = skewSymmetricMatrix(obj.data);
    end % skew
    
    function c = diff(a, b, dt)
    %DIFF Differentiate vectors.
    
      % Validate vector a attributes
      validateattributes(a, ...
        {'Vector3d'}, ...
        {'scalar'}, ...
        'Vector3d.diff', ...
        'a');
      
      % Validate vector b attributes
      validateattributes(b, ...
        {'Vector3d'}, ...
        {'scalar'}, ...
        'Vector3d.diff', ...
        'b');
      
      % Validate time step attributes
      validateattributes(dt, ...
        {'double'}, ...
        {'scalar'}, ...
        'Vector3d.diff', ...
        'dt');
      
      % Differentiate
      c = (b - a)/dt;
    end % diff
  end % methods

  % GET/SET METHODS =======================================================
  methods
    % Vector --------------------------------------------------------------
    function value = getValue(obj)
    %GETVALUE Get value.
      value = obj.data;
    end % getValue
    
    function x = getX(obj)
    %GETX Get vector x-coordinate.
      x = obj.data(1);
    end % getX
    
    function obj = setX(obj, x)
    %SETX Set vector x-coordinate.
      
      % Validate attributes
      validateattributes(x, ...
        {'double', 'sym'}, ...
        {'scalar'}, ...
        'Vector3d.setX', ...
        'x');
      
      % Set property
      obj.data(1) = x;
    end % setX

    function y = getY(obj)
    %GETY Get vector y-coordinate.
      y = obj.data(2);
    end % getY

    function obj = setY(obj, y)
    %SETY Set vector y-coordinate.
      
      % Validate attributes
      validateattributes(y, ...
        {'double', 'sym'}, ...
        {'scalar'}, ...
        'Vector3d.setY', ...
        'y');
      
      % Set property
      obj.data(2) = y;
    end % setY
    
    function z = getZ(obj)
    %GETZ Get vector z-coordinate.
      z = obj.data(3);
    end % getZ

    function obj = setZ(obj, z)
    %SETZ Set vector z-coordinate.
      
      % Validate attributes
      validateattributes(z, ...
        {'double', 'sym'}, ...
        {'scalar'}, ...
        'Vector3d.setZ', ...
        'z');
      
      % Set property
      obj.data(3) = z;
    end % setZ
  end % methods

  % STATIC METHODS ========================================================
  methods (Static = true)
    function obj = rand
    %RAND Construct uniformly distributed pseudorandom vector.
      obj = Vector3d(rand(3,1));
    end % rand
    
    function obj = randn
    %RANDN Construct normally distributed pseudorandom vector.
      obj = Vector3d(randn(3,1));
    end % randn
  end % methods
end % classdef
