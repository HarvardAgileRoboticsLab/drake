function [U] = getPotentialEnergy(obj,q)

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

U = g*m2*(l2z + c2y*cos(q2) + c2x*sin(q2)) + g*m7*(l2z - sin(q2)*(sin(q3)*(l5z - cos(q5)*(l7z + c7y*cos(q7) + c7x*sin(q7)) - sin(q5)*(c7z*sin(q6) - l7x*cos(q6) + l7y*sin(q6) + cos(q6)*(c7x*cos(q7) - c7y*sin(q7))) + l6y*cos(q5) + l6x*sin(q5)) - cos(q3)*(cos(q4)*(sin(q5)*(l7z + c7y*cos(q7) + c7x*sin(q7)) - cos(q5)*(c7z*sin(q6) - l7x*cos(q6) + l7y*sin(q6) + cos(q6)*(c7x*cos(q7) - c7y*sin(q7))) + l6x*cos(q5) - l6y*sin(q5)) - l5x*cos(q4) + l5y*sin(q4) + sin(q4)*(l6z + c7z*cos(q6) + l7y*cos(q6) + l7x*sin(q6) - sin(q6)*(c7x*cos(q7) - c7y*sin(q7)))) + l4x*cos(q3) - l4y*sin(q3)) + l3y*cos(q2) + l3x*sin(q2) + cos(q2)*(l4z - sin(q4)*(sin(q5)*(l7z + c7y*cos(q7) + c7x*sin(q7)) - cos(q5)*(c7z*sin(q6) - l7x*cos(q6) + l7y*sin(q6) + cos(q6)*(c7x*cos(q7) - c7y*sin(q7))) + l6x*cos(q5) - l6y*sin(q5)) + l5y*cos(q4) + l5x*sin(q4) + cos(q4)*(l6z + c7z*cos(q6) + l7y*cos(q6) + l7x*sin(q6) - sin(q6)*(c7x*cos(q7) - c7y*sin(q7))))) + c1z*g*m1 + g*m6*(l2z + sin(q2)*(cos(q3)*(sin(q4)*(l6z + c6y*cos(q6) + c6x*sin(q6)) + cos(q4)*(l6x*cos(q5) + c6z*sin(q5) - l6y*sin(q5) + cos(q5)*(c6x*cos(q6) - c6y*sin(q6))) - l5x*cos(q4) + l5y*sin(q4)) - l4x*cos(q3) + l4y*sin(q3) - sin(q3)*(l5z - c6z*cos(q5) + l6y*cos(q5) + l6x*sin(q5) + sin(q5)*(c6x*cos(q6) - c6y*sin(q6)))) + l3y*cos(q2) + l3x*sin(q2) + cos(q2)*(l4z + cos(q4)*(l6z + c6y*cos(q6) + c6x*sin(q6)) - sin(q4)*(l6x*cos(q5) + c6z*sin(q5) - l6y*sin(q5) + cos(q5)*(c6x*cos(q6) - c6y*sin(q6))) + l5y*cos(q4) + l5x*sin(q4))) + g*m5*(l2z - sin(q2)*(sin(q3)*(l5z + c5y*cos(q5) + c5x*sin(q5)) - cos(q3)*(c5z*sin(q4) - l5x*cos(q4) + l5y*sin(q4) + cos(q4)*(c5x*cos(q5) - c5y*sin(q5))) + l4x*cos(q3) - l4y*sin(q3)) + l3y*cos(q2) + l3x*sin(q2) + cos(q2)*(l4z + c5z*cos(q4) + l5y*cos(q4) + l5x*sin(q4) - sin(q4)*(c5x*cos(q5) - c5y*sin(q5)))) + g*m3*(l2z + c3z*cos(q2) + l3y*cos(q2) + l3x*sin(q2) - sin(q2)*(c3x*cos(q3) - c3y*sin(q3))) + g*m4*(l2z + cos(q2)*(l4z + c4y*cos(q4) + c4x*sin(q4)) - sin(q2)*(l4x*cos(q3) + c4z*sin(q3) - l4y*sin(q3) + cos(q3)*(c4x*cos(q4) - c4y*sin(q4))) + l3y*cos(q2) + l3x*sin(q2));
end