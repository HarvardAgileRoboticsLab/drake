function [G] = getKukaArmGravityVector(obj,q, c, s)

% link mass
m1=obj.m1; m2=obj.m2; m3=obj.m3; m4=obj.m4; m5=obj.m5; m6=obj.m6; m7=obj.m7;

% link 3D position
l1x=obj.l1x; l1y=obj.l1y; l1z=obj.l1z;
l2x=obj.l2x; l2y=obj.l2y; l2z=obj.l2z;
l3x=obj.l3x; l3y=obj.l3y; l3z=obj.l3z;
l4x=obj.l4x; l4y=obj.l4y; l4z=obj.l4z;
l5x=obj.l5x; l5y=obj.l5y; l5z=obj.l5z;
l6x=obj.l6x; l6y=obj.l6y; l6z=obj.l6z;
l7x=obj.l7x; l7y=obj.l7y; l7z=obj.l7z;

% gravity constant
g=obj.g;

% CoM 3D position (with repsect to joint local coordinate)
c1x = obj.c1x; c1y = obj.c1y; c1z = obj.c1z;
c2x = obj.c2x; c2y = obj.c2y; c2z = obj.c2z;
c3x = obj.c3x; c3y = obj.c3y; c3z = obj.c3z;
c4x = obj.c4x; c4y = obj.c4y; c4z = obj.c4z;
c5x = obj.c5x; c5y = obj.c5y; c5z = obj.c5z;
c6x = obj.c6x; c6y = obj.c6y; c6z = obj.c6z;
c7x = obj.c7x; c7y = obj.c7y; c7z = obj.c7z;

% joint angle data
q1 = q(1);
q2 = q(2);
q3 = q(3);
q4 = q(4);
q5 = q(5);
q6 = q(6);
q7 = q(7);

G1 = 0;
G2 = g*m2*(c2x*c(2) - c2y*s(2)) - g*m5*(c(2)*(s(3)*(l5z + c5y*c(5) + c5x*s(5)) - c(3)*(c5z*s(4) - l5x*c(4) + l5y*s(4) + c(4)*(c5x*c(5) - c5y*s(5))) + l4x*c(3) - l4y*s(3)) - l3x*c(2) + l3y*s(2) + s(2)*(l4z + c5z*c(4) + l5y*c(4) + l5x*s(4) - s(4)*(c5x*c(5) - c5y*s(5)))) - g*m3*(c3z*s(2) - l3x*c(2) + l3y*s(2) + c(2)*(c3x*c(3) - c3y*s(3))) - g*m4*(s(2)*(l4z + c4y*c(4) + c4x*s(4)) + c(2)*(l4x*c(3) + c4z*s(3) - l4y*s(3) + c(3)*(c4x*c(4) - c4y*s(4))) - l3x*c(2) + l3y*s(2)) - g*m6*(s(2)*(l4z + c(4)*(l6z + c6y*c(6) + c6x*s(6)) - s(4)*(l6x*c(5) + c6z*s(5) - l6y*s(5) + c(5)*(c6x*c(6) - c6y*s(6))) + l5y*c(4) + l5x*s(4)) - c(2)*(c(3)*(s(4)*(l6z + c6y*c(6) + c6x*s(6)) + c(4)*(l6x*c(5) + c6z*s(5) - l6y*s(5) + c(5)*(c6x*c(6) - c6y*s(6))) - l5x*c(4) + l5y*s(4)) - l4x*c(3) + l4y*s(3) - s(3)*(l5z - c6z*c(5) + l6y*c(5) + l6x*s(5) + s(5)*(c6x*c(6) - c6y*s(6)))) - l3x*c(2) + l3y*s(2)) - g*m7*(c(2)*(s(3)*(l5z - c(5)*(l7z + c7y*c(7) + c7x*s(7)) - s(5)*(c7z*s(6) - l7x*c(6) + l7y*s(6) + c(6)*(c7x*c(7) - c7y*s(7))) + l6y*c(5) + l6x*s(5)) - c(3)*(c(4)*(s(5)*(l7z + c7y*c(7) + c7x*s(7)) - c(5)*(c7z*s(6) - l7x*c(6) + l7y*s(6) + c(6)*(c7x*c(7) - c7y*s(7))) + l6x*c(5) - l6y*s(5)) - l5x*c(4) + l5y*s(4) + s(4)*(l6z + c7z*c(6) + l7y*c(6) + l7x*s(6) - s(6)*(c7x*c(7) - c7y*s(7)))) + l4x*c(3) - l4y*s(3)) - l3x*c(2) + l3y*s(2) + s(2)*(l4z - s(4)*(s(5)*(l7z + c7y*c(7) + c7x*s(7)) - c(5)*(c7z*s(6) - l7x*c(6) + l7y*s(6) + c(6)*(c7x*c(7) - c7y*s(7))) + l6x*c(5) - l6y*s(5)) + l5y*c(4) + l5x*s(4) + c(4)*(l6z + c7z*c(6) + l7y*c(6) + l7x*s(6) - s(6)*(c7x*c(7) - c7y*s(7)))));
G3 = g*m4*s(2)*(l4y*c(3) - c4z*c(3) + l4x*s(3) + s(3)*(c4x*c(4) - c4y*s(4))) - g*m6*s(2)*(s(3)*(s(4)*(l6z + c6y*c(6) + c6x*s(6)) + c(4)*(l6x*c(5) + c6z*s(5) - l6y*s(5) + c(5)*(c6x*c(6) - c6y*s(6))) - l5x*c(4) + l5y*s(4)) - l4y*c(3) - l4x*s(3) + c(3)*(l5z - c6z*c(5) + l6y*c(5) + l6x*s(5) + s(5)*(c6x*c(6) - c6y*s(6)))) - g*m7*s(2)*(s(3)*(c(4)*(s(5)*(l7z + c7y*c(7) + c7x*s(7)) - c(5)*(c7z*s(6) - l7x*c(6) + l7y*s(6) + c(6)*(c7x*c(7) - c7y*s(7))) + l6x*c(5) - l6y*s(5)) - l5x*c(4) + l5y*s(4) + s(4)*(l6z + c7z*c(6) + l7y*c(6) + l7x*s(6) - s(6)*(c7x*c(7) - c7y*s(7)))) - l4y*c(3) - l4x*s(3) + c(3)*(l5z - c(5)*(l7z + c7y*c(7) + c7x*s(7)) - s(5)*(c7z*s(6) - l7x*c(6) + l7y*s(6) + c(6)*(c7x*c(7) - c7y*s(7))) + l6y*c(5) + l6x*s(5))) - g*m5*s(2)*(c(3)*(l5z + c5y*c(5) + c5x*s(5)) + s(3)*(c5z*s(4) - l5x*c(4) + l5y*s(4) + c(4)*(c5x*c(5) - c5y*s(5))) - l4y*c(3) - l4x*s(3)) + g*m3*s(2)*(c3y*c(3) + c3x*s(3));
G4 = g*m4*(c(2)*(c4x*c(4) - c4y*s(4)) + c(3)*s(2)*(c4y*c(4) + c4x*s(4))) - g*m7*(c(2)*(c(4)*(s(5)*(l7z + c7y*c(7) + c7x*s(7)) - c(5)*(c7z*s(6) - l7x*c(6) + l7y*s(6) + c(6)*(c7x*c(7) - c7y*s(7))) + l6x*c(5) - l6y*s(5)) - l5x*c(4) + l5y*s(4) + s(4)*(l6z + c7z*c(6) + l7y*c(6) + l7x*s(6) - s(6)*(c7x*c(7) - c7y*s(7)))) - c(3)*s(2)*(l5y*c(4) - s(4)*(s(5)*(l7z + c7y*c(7) + c7x*s(7)) - c(5)*(c7z*s(6) - l7x*c(6) + l7y*s(6) + c(6)*(c7x*c(7) - c7y*s(7))) + l6x*c(5) - l6y*s(5)) + l5x*s(4) + c(4)*(l6z + c7z*c(6) + l7y*c(6) + l7x*s(6) - s(6)*(c7x*c(7) - c7y*s(7))))) - g*m5*(c(2)*(c5z*s(4) - l5x*c(4) + l5y*s(4) + c(4)*(c5x*c(5) - c5y*s(5))) - c(3)*s(2)*(c5z*c(4) + l5y*c(4) + l5x*s(4) - s(4)*(c5x*c(5) - c5y*s(5)))) - g*m6*(c(2)*(s(4)*(l6z + c6y*c(6) + c6x*s(6)) + c(4)*(l6x*c(5) + c6z*s(5) - l6y*s(5) + c(5)*(c6x*c(6) - c6y*s(6))) - l5x*c(4) + l5y*s(4)) - c(3)*s(2)*(c(4)*(l6z + c6y*c(6) + c6x*s(6)) - s(4)*(l6x*c(5) + c6z*s(5) - l6y*s(5) + c(5)*(c6x*c(6) - c6y*s(6))) + l5y*c(4) + l5x*s(4)));
G5 = -g*m6*(s(2)*(s(3)*(l6x*c(5) + c6z*s(5) - l6y*s(5) + c(5)*(c6x*c(6) - c6y*s(6))) + c(3)*c(4)*(l6y*c(5) - c6z*c(5) + l6x*s(5) + s(5)*(c6x*c(6) - c6y*s(6)))) - c(2)*s(4)*(l6y*c(5) - c6z*c(5) + l6x*s(5) + s(5)*(c6x*c(6) - c6y*s(6)))) - g*m5*(s(2)*(s(3)*(c5x*c(5) - c5y*s(5)) + c(3)*c(4)*(c5y*c(5) + c5x*s(5))) - c(2)*s(4)*(c5y*c(5) + c5x*s(5))) - g*m7*(s(2)*(s(3)*(s(5)*(l7z + c7y*c(7) + c7x*s(7)) - c(5)*(c7z*s(6) - l7x*c(6) + l7y*s(6) + c(6)*(c7x*c(7) - c7y*s(7))) + l6x*c(5) - l6y*s(5)) - c(3)*c(4)*(c(5)*(l7z + c7y*c(7) + c7x*s(7)) + s(5)*(c7z*s(6) - l7x*c(6) + l7y*s(6) + c(6)*(c7x*c(7) - c7y*s(7))) - l6y*c(5) - l6x*s(5))) + c(2)*s(4)*(c(5)*(l7z + c7y*c(7) + c7x*s(7)) + s(5)*(c7z*s(6) - l7x*c(6) + l7y*s(6) + c(6)*(c7x*c(7) - c7y*s(7))) - l6y*c(5) - l6x*s(5)));
G6 = g*m6*(c(2)*(c(4)*(c6x*c(6) - c6y*s(6)) + c(5)*s(4)*(c6y*c(6) + c6x*s(6))) + s(2)*(c(3)*(s(4)*(c6x*c(6) - c6y*s(6)) - c(4)*c(5)*(c6y*c(6) + c6x*s(6))) + s(3)*s(5)*(c6y*c(6) + c6x*s(6)))) - g*m7*(s(2)*(c(3)*(s(4)*(c7z*s(6) - l7x*c(6) + l7y*s(6) + c(6)*(c7x*c(7) - c7y*s(7))) + c(4)*c(5)*(c7z*c(6) + l7y*c(6) + l7x*s(6) - s(6)*(c7x*c(7) - c7y*s(7)))) - s(3)*s(5)*(c7z*c(6) + l7y*c(6) + l7x*s(6) - s(6)*(c7x*c(7) - c7y*s(7)))) + c(2)*(c(4)*(c7z*s(6) - l7x*c(6) + l7y*s(6) + c(6)*(c7x*c(7) - c7y*s(7))) - c(5)*s(4)*(c7z*c(6) + l7y*c(6) + l7x*s(6) - s(6)*(c7x*c(7) - c7y*s(7)))));
G7 = g*m7*(s(2)*(s(3)*(c(5)*(c7x*c(7) - c7y*s(7)) - c(6)*s(5)*(c7y*c(7) + c7x*s(7))) + c(3)*(c(4)*(s(5)*(c7x*c(7) - c7y*s(7)) + c(5)*c(6)*(c7y*c(7) + c7x*s(7))) + s(4)*s(6)*(c7y*c(7) + c7x*s(7)))) - c(2)*(s(4)*(s(5)*(c7x*c(7) - c7y*s(7)) + c(5)*c(6)*(c7y*c(7) + c7x*s(7))) - c(4)*s(6)*(c7y*c(7) + c7x*s(7))));
G = [G1;G2;G3;G4;G5;G6;G7];
end
