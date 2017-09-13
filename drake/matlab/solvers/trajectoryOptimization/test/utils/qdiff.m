function vm = qdiff(params,q1,q2,h)
vm = (q2-q1)/h;
vm(params.angle_inds) = angleDiff(q1(params.angle_inds),q2(params.angle_inds))/h;
end