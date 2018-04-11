function [f, df] = ipzt_fun(obj, act, t, x, u)

xin = [t; x; u];
[f,df] = ipzt(obj,act,xin);

% df_fd = zeros(size(df));
% step = 1e-6; %sqrt(eps(max(abs(xin))))
% dxin = step*eye(length(xin));
% for k = 1:length(xin)
%     df_fd(:,k) = (ipzt(obj,act,xin+dxin(:,k)) - ...
%         ipzt(obj,act,xin-dxin(:,k)))/(2*step);
% end


% disp('Inverse pzt model error:');
% disp(max(abs(df_fd(:)-df(:))));
% max(abs(df(:)))


end

function [V, dV] = ipzt(obj, act, xin)

nq = obj.getNumActuatedDOF();
nu = nq; %obj.getNumInputs();

t = xin(1);
qact = xin(1+(1:nq));
u = xin(1+nq+(1:nu));

V = zeros(numel(u), 1);
dVdx = zeros(numel(u), numel(qact));
dVdu = zeros(numel(u));

for i = 1:numel(qact)
    obj = act.dummy_bender(i); 
    Ct = (0.75/obj.ap.lact)*(obj.ap.f31_eff_b + (-obj.orien*qact(i)/obj.df)*...
        (obj.ap.f31_eff_f - obj.ap.f31_eff_b))*obj.ap.wn*(obj.ap.tpzt+obj.tcf)*obj.GF;
    Cb = (0.75/obj.ap.lact)*(obj.ap.f31_eff_b + (obj.orien*qact(i)/obj.df)*...
        (obj.ap.f31_eff_f - obj.ap.f31_eff_b))*obj.ap.wn*(obj.ap.tpzt+obj.tcf)*obj.GF;
    
    V(i) = (Ct*obj.dp.Vb + Cb*obj.dp.Vg - obj.orien*u(i))/(Ct+Cb);
    
    
    dCt_dq = (0.75/obj.ap.lact)*((-obj.orien/obj.df)*(obj.ap.f31_eff_f - obj.ap.f31_eff_b))*obj.ap.wn*(obj.ap.tpzt+obj.tcf)*obj.GF;
    dCb_dq = (0.75/obj.ap.lact)*((obj.orien/obj.df)*(obj.ap.f31_eff_f - obj.ap.f31_eff_b))*obj.ap.wn*(obj.ap.tpzt+obj.tcf)*obj.GF;
    
    d1dq = obj.dp.Vb*dCt_dq/(Ct + Cb);
    d2dq = obj.dp.Vg*dCb_dq/(Ct + Cb);
    d3dq = (u(i)*obj.orien*(dCt_dq + dCb_dq))/(Ct + Cb)^2;

    dVdx(i,i) = d1dq + d2dq + d3dq;    
    
    %d/du
    dVdu(i,i) = -obj.orien/(Ct + Cb);
end

dV = [zeros(numel(V), 1), dVdx, dVdu];
end
