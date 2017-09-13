function [dT,d2T] = kineticEnergyDerivatives(obj, q, v)
% kineticEnergyDerivatives  Calculate gradient and Hessian of kinetic energy.
% [dT,d2T] = kineticEnergyDerivatives(obj, q, v)
% calculates the gradient and Hessian of the system's kinetic energy with
% respect to [q; v].

nq = length(q);
nv = length(v);

options.use_mex = true;
options.compute_gradients = true;
kinsol_v = doKinematics(obj, q, v, options);

T = kineticEnergyGradmex(obj.mex_model_ptr, kinsol_v.mex_ptr);
[dT, d2T] = eval(T);

end

