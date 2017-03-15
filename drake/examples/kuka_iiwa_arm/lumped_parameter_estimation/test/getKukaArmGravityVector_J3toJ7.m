function [G] = getKukaArmGravityVector_J3toJ7(obj,q, c, s)

% link mass
m3=obj.m3; m4=obj.m4; m5=obj.m5; m6=obj.m6; m7=obj.m7;

% link 3D position
l3x=obj.l3x; l3y=obj.l3y; l3z=obj.l3z;
l4x=obj.l4x; l4y=obj.l4y; l4z=obj.l4z;
l5x=obj.l5x; l5y=obj.l5y; l5z=obj.l5z;
l6x=obj.l6x; l6y=obj.l6y; l6z=obj.l6z;
l7x=obj.l7x; l7y=obj.l7y; l7z=obj.l7z;

% gravity constant
g=obj.g;

% CoM 3D position (with repsect to joint local coordinate)
c3x = obj.c3x; c3y = obj.c3y; c3z = obj.c3z;
c4x = obj.c4x; c4y = obj.c4y; c4z = obj.c4z;
c5x = obj.c5x; c5y = obj.c5y; c5z = obj.c5z;
c6x = obj.c6x; c6y = obj.c6y; c6z = obj.c6z;
c7x = obj.c7x; c7y = obj.c7y; c7z = obj.c7z;

% joint angle data

q3 = q(1);
q4 = q(2);
q5 = q(3);
q6 = q(4);
q7 = q(5);

G1 = 0;
G2 = - (67*g*m4*s(2))/1000 - g*m6*((369*s(2))/2000 + s(2)*((3*c(4))/5000 + 431/2000) + c(2)*(s(3)/2500 - (3*c(3)*s(4))/5000)) - g*m7*((2*s(2))/5 + (101*c(4)*s(2))/1000 - (101*c(2)*c(3)*s(4))/1000) - g*m5*((521*s(2))/2000 + c(2)*(c(3)/10000 - (21*s(3))/1000));
G3 = g*m5*s(2)*((21*c(3))/1000 + s(3)/10000) - g*m6*s(2)*(c(3)/2500 + (3*s(3)*s(4))/5000) - (101*g*m7*s(2)*s(3)*s(4))/1000;
G4 = - g*m7*((101*c(2)*s(4))/1000 - (101*c(3)*c(4)*s(2))/1000) - g*m6*((3*c(2)*s(4))/5000 - (3*c(3)*c(4)*s(2))/5000);
G5 = 0;
G = [G1;G2;G3;G4;G5];
end