%CASSIESIMULATIONPARAMETERS Cassie simulation parameters.
%
% TODO:
%   - Add external force functions in here
%   - Add method for setting intial conditions (ik solver)
%   - Add flag for data logging?

% Copyright 2016-2017 Mikhail S. Jones

classdef CassieSimulationParameters %#codegen

  % CONSTANT PROPERTIES ===================================================
  properties (Constant)
    % GENERAL -------------------------------------------------------------
    
    % Acceleration of gravity [x,y,z] (m/s^2)
    gravity = [0; 0; -9.81]
    % Start time (sec)
    startTime = 0
    % Stop time (sec)
    stopTime = 60
    
    
    % VISUAL --------------------------------------------------------------
    
    % True color of robot under direct white light [r,g,b,a] (0-1)
    diffuseColor = [1 1 1 1]
    % Color of specular highlights [r,g,b,a] (0-1)
    specularColor = [0.5 0.5 0.5 1]
    % Color of shadow areas in diffuse ambient light [r,g,b,a] (0-1)
    ambientColor = [0.15 0.15 0.15 1]
    % Surface color due to self illumination [r,g,b,a] (0-1)
    emissiveColor = [0 0 0 1]
    % Sharpness of specular light reflections (0-128)
    shininess = 75
  end % properties
end % classdef