function [xA, xB, idxA, idxB] = getCollisions(obj,kinsol)
    persistent bodyAPoints bodyBPoints bodyIndexA bodyIndexB
    
    if isempty(bodyAPoints) || ~obj.isPointContact
        [~, ~, bodyAPoints, bodyBPoints, bodyIndexA, bodyIndexB] = obj.collisionDetect(kinsol);
    end
    
    xA = bodyAPoints;
    xB = bodyBPoints;
    idxA = bodyIndexA;
    idxB = bodyIndexB;
end