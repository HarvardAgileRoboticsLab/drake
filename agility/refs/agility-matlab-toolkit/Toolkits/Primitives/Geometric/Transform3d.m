%TRANSFORM3D Represents homogeneous transform in three dimensional space.
%
% Syntax:
%   transform = Transform3d();
%   transform = Transform3d(eye(4));
%   transform = Transform3d(Vector3d);
%   transform = Transform3d(Rotation3d);
%   transform = Transform3d(Vector3d, Rotation3d);
%
% Notes:
%   - Leverages orthonormal property of rotation matrices to reduce
%   computation time (inv(T) == T')

% Copyright 2016-2017 Mikhail S. Jones

classdef Transform3d %#codegen

  % PROTECTED PROPERTIES ==================================================
  properties (Access = protected)
    % Position vector object
    position
    % Rotation matrix object
    rotation
  end % properties

  % PUBLIC METHODS ========================================================
  methods (Access = public)
    function obj = Transform3d(varargin)
    %TRANSFORM3D Construct a transformation in three dimensions.

      % Set code generation settings
      coder.allowpcode('plain');
      
      % Check the number of input arguments
      switch nargin
        case 0
          % Initialize position vector and rotation matrix
          obj.position = Vector3d;
          obj.rotation = Rotation3d;

        case 1
          % Check input argument type
          switch class(varargin{1})
            case 'Transform3d'
              % Already a transform object
              obj = varargin{1};
            case 'Vector3d'
              % Set user specified position vector
              obj.position = varargin{1};
              obj.rotation = Rotation3d;

            case 'Rotation3d'
              % Set user specified rotation matrix
              obj.position = Vector3d;
              obj.rotation = varargin{1};

            otherwise
              % Validate input argument attributes
              validateattributes(varargin{1}, ...
                {'double', 'sym'}, ...
                {'real', 'size', [4, 4]}, ...
                'Transform3d', ...
                'transformationMatrix');

              % Set user specified matrix
              obj.position = Vector3d(varargin{1}(1:3,4));
              obj.rotation = Rotation3d(varargin{1}(1:3,1:3));
          end % switch

        case 2
          % Validate position attributes
          validateattributes(varargin{1}, ...
            {'Vector3d'}, ...
            {'scalar'}, ...
            'Transform3d', ...
            'position');

          % Validate rotation attributes
          validateattributes(varargin{2}, ...
            {'Rotation3d'}, ...
            {'scalar'}, ...
            'Transform3d', ...
            'rotation');

          % Set user specified vectors
          obj.position = varargin{1};
          obj.rotation = varargin{2};

        otherwise
          % Throw error
          error('Too many input arguments.');
      end % switch
    end % Transform3d
    
    function plot(obj)
    %PLOT Plot the coordinate system for the transform.
    
      % Scale of frame representation
      scale = 0.1;
      
      % X-axis translation
      tx = plot3([0 scale], [0 0], [0 0], ...
        'Color', 'r', ...
        'LineWidth', 3);
      
      % Y-axis translation
      ty = plot3([0 0], [0 scale], [0 0], ...
        'Color', 'g', ...
        'LineWidth', 3);

      % Z-axis translation
      tz = plot3([0 0], [0 0], [0 scale], ...
        'Color','b', ...
        'LineWidth', 3);
      
      % Construct handle graphics transform object
      hg = hgtransform;

      % Set plots to be children of handle graphics transform
      set(tx, 'Parent', hg);
      set(ty, 'Parent', hg);
      set(tz, 'Parent', hg);

      % Set handle graphics object transformation matrix
      set(hg, 'Matrix', obj.getValue);
    end % plot
    
    function disp(obj)
    %DISP Display object value.
      disp(obj.getValue);
    end % disp
    
    % MATH OPERATOR METHODS -----------------------------------------------
    function c = mtimes(a, b)
    %MTIMES Matrix multiplication.

      % Validate transform A attributes
      validateattributes(a, ...
        {'Transform3d'}, ...
        {'scalar'}, ...
        'Transform3d.mtimes', ...
        'a');
      
      % Validate transform B attributes
      validateattributes(b, ...
        {'Transform3d'}, ...
        {'scalar'}, ...
        'Transform3d.mtimes', ...
        'b');
      
      % Compute position vector
      p = a.position + a.rotation*b.position;

      % Compute rotation matrix
      R = a.rotation*b.rotation;

      % Construct transform
      c = Transform3d(p, R);
    end % mtimes

    function c = mrdivide(a, b)
    %MRDIVIDE Matrix right handed division.

      % Validate transform A attributes
      validateattributes(a, ...
        {'Transform3d'}, ...
        {'scalar'}, ...
        'Transform3d.mrdivide', ...
        'a');
      
      % Validate transform B attributes
      validateattributes(b, ...
        {'Transform3d'}, ...
        {'scalar'}, ...
        'Transform3d.mrdivide', ...
        'b');
      
      % Compute position vector
      p = a.position - a.rotation*(b.rotation\b.position);

      % Compute rotation matrix
      R = a.rotation/b.rotation;

      % Construct transform
      c = Transform3d(p, R);
    end % mrdivide

    function c = mldivide(a, b)
    %MLDIVIDE Matrix left handed division.

      % Validate transform A attributes
      validateattributes(a, ...
        {'Transform3d'}, ...
        {'scalar'}, ...
        'Transform3d.mldivide', ...
        'a');
      
      % Validate transform B attributes
      validateattributes(b, ...
        {'Transform3d'}, ...
        {'scalar'}, ...
        'Transform3d.mldivide', ...
        'b');
      
      % Compute rotation matrix
      R = a.rotation\b.rotation;

      % Compute position vector
      p = a.rotation\b.position - a.rotation\a.position;

      % Construct transform
      c = Transform3d(p, R);
    end % mldivide
    
    function obj = inverse(obj)
    %INVERSE Matrix inverse.
    
      % Get local variables
      R_inv = obj.rotation.inverse;
      p = obj.position;
      
      % Compute rotation matrix inverse
      obj.rotation = R_inv;
      
      % Compute position vector expressed in new rotation
      obj.position = -R_inv*p;
    end % inverse
    
    % TRANSFORMATION OPERATOR METHODS -------------------------------------
    function obj = translate(obj, xyz)
    %TRANSLATE Translate transform.
      obj.position = obj.position + obj.rotation*Vector3d(xyz);
    end % translate
    
    function obj = rotX(obj, varargin)
    %ROTX Rotate around X-axis.
      obj.rotation = obj.rotation.rotX(varargin{:});
    end % rotX

    function obj = rotY(obj, varargin)
    %ROTY Rotate around Y-axis.
      obj.rotation = obj.rotation.rotY(varargin{:});
    end % rotY

    function obj = rotZ(obj, varargin)
    %ROTZ Rotate around Z-axis.
      obj.rotation = obj.rotation.rotZ(varargin{:});
    end % rotZ

    function obj = rotZYX(obj, varargin)
    %ROTZYX Rotate using ZYX Euler angles.
      obj.rotation = obj.rotation.rotZYX(varargin{:});
    end % rotZYX

    function obj = rotZYZ(obj, varargin)
    %ROTZYZ Rotate using ZYZ Euler angles.
      obj.rotation = obj.rotation.rotZYZ(varargin{:});
    end % rotZYZ
    
    function obj = rotAxis(obj, varargin)
    %ROTAXIS Rotate about arbitrary axis.
      obj.rotation = obj.rotation.rotAxis(varargin{:});
    end % rotAxis
  end % methods

  % GET/SET METHODS =======================================================
  methods
    % Transformation ------------------------------------------------------
    function value = getValue(obj)
    %GETVALUE Get value.
      value = [obj.rotation.getValue, obj.position.getValue; 0, 0, 0, 1];
    end % getValue
    
    function position = getPosition(obj)
    %GETPOSITION Get position vector from transform.
      position = obj.position;
    end % getPosition

    function rotation = getRotation(obj)
    %GETROTATION Get rotation matrix from transform.
      rotation = obj.rotation;
    end % getRotation
    
    % Adjoint Transformations ---------------------------------------------
    function rotationAdjoint = getRotationAdjoint(obj)
    %GETROTATIONADJOINT Get the rotation adjoint transformation matrix.
    %
    % Notes:
    %   V' = AdR*V
    %
    %   v' = R*v
    %   w' = R*w
    
      % Get local variables 
      R = obj.rotation.getValue;
      
      % Construct the rotation adjoint transformation matrix
      rotationAdjoint = [...
        R, zeros(3); ...
        zeros(3), R];
    end % getRotationAdjoint
    
    function adjoint = getAdjoint(obj)
    %GETADJOINT Get the adjoint transformation matrix.
    %
    % Notes:
    %   V' = AdT*V
    %
    %   v' = R*v + p x R*w
    %   w' = R*w
    
      % Get local variables
      R = obj.rotation.getValue;
      p_x = obj.position.skew;
      
      % Construct the adjoint transformation matrix
      adjoint = [...
        R, p_x*R; ...
        zeros(3), R];
    end % getAdjoint
    
    function adjointInverse = getAdjointInverse(obj)
    %GETADJOINTINVERSE Get the adjoint inverse transformation matrix.
    
      % Get local variables
      R_inv = obj.rotation.inverse.getValue;
      p_x = obj.position.skew;
      
      % Construct the adjoint inverse transformation matrix
      adjointInverse = [...
        R_inv, -R_inv*p_x; ...
        zeros(3), R_inv];
    end % getAdjointInverse
    
    function adjointTranspose = getAdjointTranspose(obj)
    %GETADJOINTTRANSPOSE Get the transposed adjoint transformation matrix.
    
      % Get local variables
      R_T = obj.rotation.transpose.getValue;
      p_x = obj.position.skew;
      
      % Construct the transposed adjoint inverse transformation matrix
      adjointTranspose = [...
        R_T, zeros(3); ...
        R_T*p_x.', R_T];
    end % getAdjointTranpose
    
    function adjointInverseTranspose = getAdjointInverseTranspose(obj)
    %GETADJOINTINVERSETRANSPOSE Get the transpose adjoint inverse transformation matrix.
    
      % Get local variables
      R = obj.rotation.getValue;
      p_x = obj.position.skew;
      
      % Construct the transposed adjoint transformation matrix
      adjointInverseTranspose = [...
        R, zeros(3); ...
        -p_x.'*R, R];
    end % getAdjointInverseTranpose
  end % methods
    
  % STATIC METHODS ========================================================
  methods (Static = true)
    function obj = rand
    %RAND Construct uniformly distributed pseudorandom transform.
      
      % Initialize vector and rotation objects
      vector = Vector3d;
      rotation = Rotation3d;
      
      % Construct random transform object
      obj = Transform3d(vector.rand, rotation.rand);
    end % rand
    
    function obj = randn
    %RANDN Construct normally distributed pseudorandom transform.
      
      % Initialize vector and rotation objects
      vector = Vector3d;
      rotation = Rotation3d;
      
      % Construct random transform object
      obj = Transform3d(vector.randn, rotation.randn);
    end % randn
  end % methods
end % classdef
