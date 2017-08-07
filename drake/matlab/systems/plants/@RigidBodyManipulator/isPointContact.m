function pt_geometry = isPointContact(obj)
    pt_geometry = true;
    
    for i=2:numel(obj.body) % iterate through all of the bodies
        body = obj.body(i);
        for j = 1:numel(body.collision_geometry) % iterate through all of the collision geometries
            if ~isa(body.collision_geometry{i},'RigidBodySphere') || body.collision_geometry{i}.radius > 1e-5
                % exit as soon as it fails
                pt_geometry = false;
                return;
            end
        end
    end
end