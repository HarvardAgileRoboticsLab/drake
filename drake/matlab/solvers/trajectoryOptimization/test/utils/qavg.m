function qm = qavg(params,q1,q2)

angle_inds = params.angle_inds; 

qm = (q1+q2)/2;
if ~isempty(angle_inds)
    qm(angle_inds) = angleAverage(q1(angle_inds),q2(angle_inds));
end
end