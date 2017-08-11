function [f, df] = forwardDyn(zin, params)

plant = params.plant;
h1 = params.h; 
h2 = params.h; 
q1 = params.q1;
q2 = params.q2; 
kl1 = params.kl1; 

nKL = plant.getNumStateConstraints()-4;  % NDD: constraint forces
nQ = plant.getNumPositions();

q3 = zin(1:nQ);
kl2 = zin(nQ+(1:nKL));


%Take care of angle wrap-around
qm1 = qavg(params,q1,q2);
vm1 = qdiff(params,q1,q2,h1);
qm2 = qavg(params,q2,q3);
vm2 = qdiff(params,q2,q3,h2);

%Discrete Euler-Lagrange equation
[D1L1,D2L1,~,~,~,~,~] = LagrangianDerivs(params,qm1,vm1);
[D1L2,D2L2,D1D1L2,D1D2L2,D2D2L2,~,~] = LagrangianDerivs(params,qm2,vm2);

f_del = (h1/2)*D1L1 + D2L1 + (h2/2)*D1L2 - D2L2;
df_del = [(h2/2)*((1/2)*D1D1L2+(1/h2)*D1D2L2') - (1/2)*D1D2L2 - (1/h2)*D2D2L2, ... % d/dq3
    zeros(nQ, nKL)]; %d/dkl

%NDD: closed chains
good_inds = [1; 3]; 
[~, dKC1, ~] = plant.positionConstraints(qm1);
[~, dKC2, dKCdqm2] = plant.positionConstraints(qm2);

dKC1 = dKC1(good_inds, :); 
dKC2 = dKC2(good_inds, :); 
dKCdqm2 = dKCdqm2(good_inds, :);
dKCdqm2 = reshape(dKCdqm2', nQ, nKL*nQ)';

%Total dynamics residual incluing control + contact forces
f_dyn = f_del + (h1/2)*dKC1'*kl1 + (h2/2)*dKC2'*kl2;

%Dynamics Derivatives
df_dyn = df_del + [(h2/4)*kron(kl2', eye(nQ))*dKCdqm2, ...% d/dq3
    (h2/2)*dKC2']; %d/dkl2

[phi, dphi, ~] = plant.positionConstraints(q3); 
phi = phi(good_inds); 
dphi = dphi(good_inds, :); 
f = [f_dyn; phi]; 
df = [df_dyn; [dphi, zeros(nKL)]]; 
end





