%UDPINITFCN UDP communication initialization callback function.
%
% Syntax:
%   udpInitFcn;

% Copyright 2017 Agility Robotics
% Author: Daniel Bennett

% Initialize Cassie input structure byte packing information
[inputBytes, inputTypes, inputSizes] = cassieInputs.getUdpPackingInfo;

% Initialize Cassie output structure byte packing information
[outputBytes, outputTypes, outputSizes] = cassieOutputs.getUdpPackingInfo;