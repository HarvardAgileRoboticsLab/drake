 function [f, df] = SLFitFun(z, params)

SL = params.SL; 
nq = SL.getNumPositions();

SLSimple = params.SLSimple; 
nqS = SLSimple.getNumPositions();
nvS = SLSimple.getNumVelocities(); 

t = params.t; 
x = params.x; 
u = params.u; 
tau = params.tau; 
sfb_inputs = params.sfb_inputs; 

k = z(1:nqS^2); 
p = z(nqS^2+(1:nvS^2)); 

f = zeros(numel(t)*nqS,1);
df = zeros(numel(t)*nqS, numel(z)); 
h = SL.timestep; 

for i = 1:numel(t) 
    
    xip1_full = SL.update(t(i), x(:,i), u(:,i));    
    
    xi_simple = [x(sfb_inputs,i) + SLSimple.q0; x(nq+sfb_inputs,i)];
    
    [xidot_simple, dK, dP] =  SLSimple.dynamics_sd_fit_fun([t(i); xi_simple; tau(:,i); k; p]);
        
    xip1_simple = xi_simple +xidot_simple*h;
    
    f(nqS*(i-1)+1:(nqS*i)) = xip1_full(nq+sfb_inputs) - xip1_simple(nqS+(1:nvS));        
    
    df(nqS*(i-1)+1:(nqS*i), :) = -h*[dK(nqS+(1:nvS),:), dP(nqS+(1:nvS), :)]; 
%     f(nq*(i-1)+1:(nq*i)) = (1/2)*erri'*erri;
%     df(i,:) = (h*xi_simple(nqS+(1:nvS)) - h*xip1_full(nq+sfb_inputs) + h^2*xidot_simple(nqS+(1:nvS)))'*...
%         [dK(nqS+(1:nvS),:), dP(nqS+(1:nvS), :)];     
end



end



