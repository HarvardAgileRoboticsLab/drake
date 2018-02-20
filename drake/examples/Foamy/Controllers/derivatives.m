function d = derivatives(t, phi)
    dt = t(2) - t(1);
    d = zeros(size(phi));
    d(2) = (-1/2 * phi(1) + 1/2 * phi(3))/dt;
    
    for i = 3 : numel(t) - 2
        d(i) = 1/dt * (1/12 * phi(i-2) - 2/3 * phi(i-1) + 2/3 * phi(i+1) - 1/12 * phi(i+2));
    end
    
    d(end - 1) = (-1/2 * phi(end-2) + 1/2 * phi(end))/dt;
    d(1) = d(2);
    d(end) = d(end-1);
end