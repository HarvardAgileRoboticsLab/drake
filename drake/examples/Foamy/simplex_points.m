function chi = simplex_points(n)

Np = n+2;
w = 1/Np;

chi = zeros(n,Np);
chi(1,2) = -1/sqrt(2*w);
chi(1,3) = 1/sqrt(2*w);
for k = 2:n
    chi(k,1+(1:k)) = -1/sqrt(k*(k+1)*w);
    chi(k,2+k) = k/sqrt(k*(k+1)*w);
end

end

