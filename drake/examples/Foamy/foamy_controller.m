function u = foamy_controller(x,x0,u0,K)

q = x(4:7);
q0 = x0(4:7);
v = [zeros(3,1), eye(3)]*qmultiply(qconj(q0),q);

u = u0 - K*[x(1:3)-x0(1:3); v; x(8:13)-x0(8:13)]; 

end

% ---------- Helper Functions ---------- %
function q = qmultiply(q1, q2)
    %Quaternion multiplication
    s1 = q1(1);
    v1 = q1(2:4);
    s2 = q2(1);
    v2 = q2(2:4);
    q = [s1*s2 - v1'*v2; s1*v2 + s2*v1 + cross(v1,v2)];
end

function qc = qconj(q)
    %Quaternion conjugate
    qc = [q(1); -q(2:4)];
end
