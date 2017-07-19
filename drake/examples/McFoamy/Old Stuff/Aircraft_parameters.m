function [PropGeom,Geom,CG] = Aircraft_parameters()
%MCFOAMY GEOMETRY

%        x     y    z
CG = [-270, 0, 5.89]; %measured from propeller plane

PropGeom = [12, 254];

%         #segs  wing root   wing tip  wing span
WingGeom = [7,    260.41,    152.40,    431.80, 0, 0, 0];
TailGeom = [3,    152.73,    139.70,    181.61, 0, 0, 0];
RudGeom =  [4,    204.01,    128.55,    215.90, 0, 0, 0];
BodyGeom = [4,    641.71,    641.71,    146.27, 0, 0, 0];

% Wing Sections 
%            Span    Area     Chord FlapChord   x     y     z
WingSec1 = [57.26, 14499.66, 253.40,      0, -217.07,  33.38, 0];
WingSec2 = [58.72, 13977.71, 237.58, 102.04, -211.39,  91.76, 0];
WingSec3 = [59.29, 13198.25, 222.69,  99.23, -215.43, 151.14, 0];
WingSec4 = [62.58, 12974.40, 207.42,  96.36, -212.18, 212.01, 0];
WingSec5 = [62.58, 11992.21, 191.74,  93.41, -208.85, 274.55, 0];
WingSec6 = [63.25, 11122.51, 175.97,  90.44, -205.50, 337.42, 0];
WingSec7 = [61.90,   9913.9, 160.28,  87.49, -202.26, 400.21, 0];

% Tail Sections 
%            Span    Area     Chord FlapChord   x     y     z
TailSec1 = [66.55, 7241.76, 117.41,  81.87, -711.27,  39.53, 0];
TailSec2 = [51.09, 7452.58, 145.64, 111.49, -719.72,  93.33, 0];
TailSec3 = [62.02, 8892.67, 143.16, 143.16, -715.15, 151.57, 0];

% Rudder Sections
%            Span    Area     Chord FlapChord   x     y     z
RudSec1 = [46.30,  9291.95, 200.62, 120.36, -733.94, 0,   23.54];
RudSec2 = [62.19, 11997.38, 192.95, 112.69, -735.81, 0,  -29.96];
RudSec3 = [69.31, 11715.82, 169.49, 103.43, -744.14, 0,  -94.63];
RudSec4 = [38.10,  5298.95, 141.14, 141.14, -760.76, 0, -148.56];

% Body Sections
%            Span    Area    Chord  FlapChord   x     y     z
BodySec1 = [35.87, 23018.14, 641.71,   0,   -205.09,  0,  53.81];
BodySec2 = [35.87, 23018.14, 641.71,   0,   -205.09,  0,  17.94];
BodySec3 = [37.26, 23910.11, 641.71,   0,   -205.09,  0, -18.63];
BodySec4 = [37.26, 23910.11, 641.71,   0,   -205.09,  0, -55.89];

Geom = [WingGeom; TailGeom; RudGeom; BodyGeom;...
        WingSec1; WingSec2; WingSec3; WingSec4; WingSec5; WingSec6; WingSec7;...
        TailSec1; TailSec2; TailSec3; ...
        RudSec1; RudSec2; RudSec3; RudSec4;...
        BodySec1; BodySec2; BodySec3; BodySec4]';
end

