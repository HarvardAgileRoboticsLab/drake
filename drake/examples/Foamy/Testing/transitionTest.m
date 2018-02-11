

function transitionTest()


    p = HybridFoamyPendulumPlant(.1);
    %x = [0;0;0;3;0;0;0;0;0;0;0;0;0];
    x = [0,0,0,-0.000311716848117638,0.990660551393713,-0.00226509797146184,0.136332109483025,0.0135059548822321,0,-0.348141345875652,0,1,0,0,6,0,0,0,0,0,6,0,0,0,0,0];
    x = x';
    u = [1;1;1;1];
    %[~,dc] = contact_Force_Phi_Constraint(k,x,u,N,p,nq);
    t = 1;
    [~,~,~,dxp] = transitionT(p,1,t,x,u);
    dxp = dxp(:,2:end);
    options.grad_method = 'numerical';
    [~,dxp_test] = geval(@(t,x,u)transitionT(p,1,t,x,u),t,x,u,options);
    valuecheck(dxp,dxp_test,1e-5)

end


            function [to_mode_xn,to_mode_num,status,dxp] = transitionT(obj,from_mode_num,t,mode_x,u)
                mode_x = double(mode_x);
                u = double(u);
                [xdot,dxdot] = impact_dynamicsT(obj,t,mode_x,u);
                dt = 0.01;
                to_mode_xn = mode_x + dt*xdot;
                %to_mode_xn = mode_x;
                to_mode_num = 2;
                status = 0;
                %[tempa,tempb] = obj.modes{from_mode_num}.dynamics(t,[mode_x],u);
                %dxp = obj.dynamics(t,[from_mode_num;mode_x],u);
                %dxp = [zeros(13,1),tempb];
                
                dxp_pre = [zeros(26,1),dxdot];
                dxp_pre = dt*dxp_pre + [zeros(26,2),eye(26),zeros(26,4)];
                
                dxp = dxp_pre;
                %dxp = [zeros(26,2),eye(26),zeros(26,4)];
            end
            
              function [xdot, dxdot] = impact_dynamicsT(obj,t,x,u)
            F_impact = 5;
            if nargout == 1
                %xdot = foamy_dynamics(t,x,u);
                xdot = impact_foamy_pendulum_dynamics_mex(t,x,u,obj.m_package,F_impact);
                dxdot = 0;
            else %Need Jacobian
                options.grad_level = 0;
                options.grad_method = 'numerical';
                options.diff_type = 'central';
                %[xdot, dxdot] = geval(@(t1,x1,u1) foamy_dynamics(t1,x1,u1),t,x,u,options);
                [xdot, dxdot] = geval(@(t1,x1,u1) impact_foamy_pendulum_dynamics_mex(t1,x1,u1,obj.m_package,F_impact),t,x,u,options);
            end
        end