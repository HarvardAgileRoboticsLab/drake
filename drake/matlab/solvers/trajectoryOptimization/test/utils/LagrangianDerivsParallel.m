function [D1L,D2L,D1D1L,D1D2L,D2D2L,B,dBdq] = LagrangianDerivsParallel(obj,q2,v)
nq = length(q2);
nv = length(v);
[M,G,B,dM,dG,dB] = manipulatorDynamics(obj, q2, zeros(nv,1));

dM = reshape(dM,nq*nq,nq+nv);
dMdq = dM(:,1:nq);
dBdq = dB(:,1:nq);

D1L = 0.5*dMdq'*kron(v,v) - G;
D2L = M*v;

%D1D1L = -dG(:,1:nq); %throwing out second derivative of M terms here

D1D1L = zeros(nq);
% deltaq = step*eye(nq);


ticBytes(gcp); 
parfor k = 1:nq
    step = sqrt(eps(max(q2)));
    deltaq = zeros(nq,1); 
    deltaq(k) = step; 
    
    [~,Gp,~,dMp] = manipulatorDynamics(obj, q2+deltaq, zeros(nv,1));
    dMp = reshape(dMp,nq*nq,nq+nv);
    dMdqp = dMp(:,1:nq);
    
    [~,Gm,~,dMm] = manipulatorDynamics(obj, q2-deltaq, zeros(nv,1));
    dMm = reshape(dMm,nq*nq,nq+nv);
    dMdqm = dMm(:,1:nq);
    
    D1p = 0.5*dMdqp'*kron(v,v) - Gp;
    D1m = 0.5*dMdqm'*kron(v,v) - Gm;
    
    D1D1L(:,k) = (D1p - D1m)/(2*step);
end
tocBytes(gcp)

%disp(sprintf('D1D1L error: %d',max(abs(D1D1L_fd(:)-D1D1L(:)))));

D1D2L = kron(v',eye(nq))*dMdq;
D2D2L = M;

end