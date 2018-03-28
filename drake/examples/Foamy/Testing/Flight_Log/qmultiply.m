function q = qmultiply(q1,q2)

    %Quaternion multiplication
    s1 = q1(1);
    v1 = q1(2:4);
    s2 = q2(1);
    v2 = q2(2:4);
    q = [s1*s2 - v1'*v2; s1*v2 + s2*v1 + cross(v1,v2)];

end


%0 = roll 180,yaw 180: pitch 180
%180x =[0 1 0 0] = [-1 0 3 -2] CLOSE or [-1 0 -3 2] NOPE
%180y =[0 0 1 0] = [-2 -3 0 1] NOPE or [-2 3 0 -1] NOPE
%180z =[0 0 0 1] = [-3 2 -1 0] NOPE or [-3 -2 1 0] NOPE