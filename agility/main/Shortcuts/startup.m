%STARTUP This script defines a project startup shortcut.

% Copyright 2017 Agility Robotics
% Author: Mikhail Jones

% Open project object for editing via command line
project = simulinkproject();

% Get project root folder
projectRoot = project.RootFolder;

% Construct paths to cache and code folders
myCacheFolder = fullfile(projectRoot, 'Build', 'Cache');
myCodeFolder = fullfile(projectRoot, 'Build', 'Code');

% Set the file generation folder paths
Simulink.fileGenControl('set',...
    'CacheFolder', myCacheFolder,...
    'CodeGenFolder', myCodeFolder,...
    'createDir', true);