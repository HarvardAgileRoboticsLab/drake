function Quat = EulToQuat(Eul)

phi = Eul(1);
the = Eul(2);
psi = Eul(3);

Quat = zeros(1,4);

Quat(1) = cos(psi/2)*cos(the/2)*cos(phi/2) + sin(psi/2)*sin(the/2)*sin(phi/2);
Quat(2) = cos(psi/2)*cos(the/2)*sin(phi/2) - sin(psi/2)*sin(the/2)*cos(phi/2);
Quat(3) = cos(psi/2)*sin(the/2)*cos(phi/2) + sin(psi/2)*cos(the/2)*sin(phi/2);
Quat(4) = sin(psi/2)*cos(the/2)*cos(phi/2) - cos(psi/2)*sin(the/2)*sin(phi/2);
