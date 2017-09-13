
function [fdamp, dfdamp] = computeDampingForcesFun(obj, xin)

% Damping at the generalized velocities
nQ = obj.getNumPositions();
v = xin;

xx = Point(obj.getStateFrame());
rb = obj.body;
for i = 1:numel(rb)
    rbi = rb(i);
    ji = obj.findPositionIndices(rbi.jointname)+nQ;
    xx(ji) = rbi.damping;
end

fdamp = -diag(xx(nQ+1:end))*v;
dfdamp = -diag(xx(nQ+1:end));


end