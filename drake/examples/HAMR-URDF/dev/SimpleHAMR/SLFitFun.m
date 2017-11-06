function [f, df] = SLFitFun(z, params)

SL = params.SL; 
nq = SL.getNumPositions();
nv = SL.getNumVelocities();
nu = SL.getNumInputs();

SLSimple = params.SLSimple; 
nqS = SLSimple.getNumPositions();
nvS = SLSimple.getNumVelocities(); 

t = params.t; 
x = params.x; 
u = params.u; 
tau = params.tau; 
sfb_inputs = params.sfb_inputs; 

k = z(1:nqS^2); 
p = z(nqs^2+(1:nqS^2)); 

f = 0;
df_dK = zeros(1, nqS*nqS); 
df_dP = zeros(1, nvS*nvS); 

for i = 1:numel(t) 
    
    xip1_full = SL.update(t(i), x(:,i), u(:,i));    
    
    xi_simple = [x(sfb_inputs,i) + SLSimple.q0; x(nq+sfb_inputs,i)];
    
    [xidot_simple, dK, dP] =  SLSimple.dynamics_sd_fit_fun([t(i); xi_simple; tau(:,i); k; p]);
    
    xip1_simple = xi_simple +xidot_simple*SL.timestep;
    
    erri = xip1_full(nqS+sfb_inputs) - xip1_simple(nqS+(1:nvS));
    
    f = f + (1/2)*erri'*erri;
    
    df_dK = df_dK + erri'*dK(nqS+(1:nvS), :)*SL_;    
    df_dP = df_dP + erri'*dP(nqS+(1:nvS), :);    
end

df = [df_dK, df_dP]; 

end



