

function packageGuard_Test()


    p = HybridFoamyPendulumPlant(.1);
    %x = [0;0;0;3;0;0;0;0;0;0;0;0;0];
    x = [0,10,0,-0.000311716848117638,0.990660551393713,-0.00226509797146184,0.136332109483025,0.0135059548822321,0,-0.348141345875652,0,1,0,0,6,0,0,0,0,0,6,0,0,0,0,0];
    u = [zeros(4,1)];
    %[~,dc] = contact_Force_Phi_Constraint(k,x,u,N,p,nq);
    t = 1;
    [g,dg] = packageGuardz(p,t,x,u);
    options.grad_method = 'numerical';
    [~,dg_test] = geval(@(t,x,u)packageGuardz(p,t,x,u),t,x,u,options);
    valuecheck(dg,dg_test,1e-5)

end


            function [g,dg] = packageGuardz(obj,t,x,u)
                %g = x(10) - 0.5;
                q = x(11:14);
                v = q(2:4);
                hat = [ 0  -v(3) v(2);
                v(3)  0  -v(1);
                 -v(2) v(1)  0];
                R_p = eye(3) + 2*hat*(hat + q(1)*eye(3));
                p_l = [0;0;-.3];
                p_end = x(8:10)'-R_p*p_l;

                package_location = [0;0;0];
                hoop_length = 0.25;

                pend_pack = p_end-package_location;

                g = norm(pend_pack)-hoop_length;

                %g = p_end(3)-0;
                %g = x(10);
                %dg = [0,zeros(1,9),1,zeros(1,20)];
                dpden = 1/(norm(pend_pack));
                dpdxpend = pend_pack;

                dqs = (3/5)*[q(3) -q(2) 0; q(4) -q(1) -2*q(2);...
                            q(1) q(4) -2*q(3); q(2) q(3) 0]';



                dg = [0,zeros(1,7),(dpden.*dpdxpend)', dpden*pend_pack'*dqs,zeros(1,16)];

                %dg = [0,-1,zeros(1,29)];
            end