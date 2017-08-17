function [ f, df ] = initializeDyn( zin, params )


plant = params.plant;
h = params.h;
q0 = params.q0;
v0 = params.v0;

nKL = 6*plant.getNumStateConstraints()-3;  % NDD: constraint forces
nQ = plant.getNumPositions();

q1 = zin(1:nQ);
kl = zin(nQ+(1:nKL));

%Discrete Euler-Lagrange equation
[M0,~,~,~] = manipulatorDynamics(plant, q0, zeros(nQ,1));
p0 = M0*v0;
qm = qavg(params,q0,q1);
vm = qdiff(params,q0,q1,h);
[D1L,D2L,D1D1L,D1D2L,D2D2L,~,~] = LagrangianDerivs(params,qm,vm);


f_del = p0 + (h/2)*D1L - D2L;

df_del = [(h/4)*D1D1L + (1/2)*D1D2L' - (1/2)*D1D2L - (1/h)*D2D2L, ... % d/dq1
    zeros(nQ, nKL)]; %d/dkl

% damping forces
[fdamp, dfdamp] = computeDampingForcesFun(plant, vm);

%NDD: closed chains
good_inds = [1;3;5]; %; 3; 5];
[~, dKC, dKCdqm] = plant.positionConstraints(qm);
% KC = KC(good_inds);
dKC = dKC(good_inds, :);
dKCdqm = dKCdqm(good_inds, :);
dKCdqm = reshape(dKCdqm', nQ, nKL*nQ)';

f_dyn = f_del + (h/2)*(dKC'*kl + fdamp);

df_dyn = df_del + [(h/4)*kron(kl', eye(nQ))*dKCdqm + (1/2)*dfdamp, ...  % d/dq1
    (h/2)*dKC']; % d/dkl

% f = f_dyn;
% df = df_dyn;
[phi, dphi, ~] = plant.positionConstraints(q1);
phi = phi(good_inds);
dphi = dphi(good_inds, :);
f = [f_dyn; phi];
df = [df_dyn; [dphi, zeros(nKL)]];

end

