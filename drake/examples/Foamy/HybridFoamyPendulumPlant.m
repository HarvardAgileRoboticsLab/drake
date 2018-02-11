classdef HybridFoamyPendulumPlant < HybridDrakeSystem
    
   properties
        m_package = 0;
    end
    
    methods
        
        function obj = HybridFoamyPendulumPlant(m_p)
            obj = obj@HybridDrakeSystem(4,26);
            %obj = obj@DrakeSystem(13,0,4,13,0,1);
            %obj = obj.setStateFrame(CoordinateFrame('FoamyState',13,'x',{'x1','x2','x3','q0','q1','q2','q3','v1','v2','v3','w1','w2','w3'}));
            %obj = obj.setInputFrame(CoordinateFrame('FoamyInput',4,'u',{'thr','ail','elev','rud'}));
            %obj = obj.setOutputFrame(CoordinateFrame('FoamyState',13,'x',{'x1','x2','x3','q0','q1','q2','q3','v1','v2','v3','w1','w2','w3'}));
            if nargin == 0
                obj.m_package = 0;
            else
                obj.m_package = m_p;
            end
            
            mode_1 = FoamyPendulumPlant(0);
            obj = setInputFrame(obj,mode_1.getInputFrame);
            obj = setOutputFrame(obj,mode_1.getOutputFrame);
            [obj,mode_1] = addMode(obj,mode_1,'pre-collision');
            mode_2 = FoamyPendulumPlant(obj.m_package);
            obj = setInputFrame(obj,mode_2.getInputFrame);
            obj = setOutputFrame(obj,mode_2.getOutputFrame);
            [obj,mode_2] = addMode(obj,mode_2,'post-collision');
            %guard1 = inline('-x(1)','obj','t','x','u');
            
            
            %obj = addTransition(obj,mode_1,@packageGuard,@transition,1,1)
            obj = addTransition(obj,mode_1,andGuards(obj,@packageGuardx,@packageGuardz),@transition,1,1)
            
            %obj = obj.setOutputFrame(obj.getStateFrame); %full state feedback
            %quaternion unit-norm constraint -- why doesn't the index input work?
            %obj = obj.addStateConstraint(QuadraticConstraint(.5,.5,blkdiag(zeros(3),eye(4),zeros(6)),zeros(13,1)));
        end
            
            function [g,dg] = packageGuardx(obj,t,x,u)
            g = -x(1);
            dg = [0,-1,zeros(1,29)];
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
                p_end = x(8:10)-R_p*p_l;

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
            
            function R = qtoR(q)
                 R = eye(3) + 2*hat(q(2:4))*(hat(q(2:4)) + q(1)*eye(3));
            end
            
            function C = hat(x)
            C = [ 0  -x(3) x(2);
            x(3)  0  -x(1);
             -x(2) x(1)  0];
            end
        
            function [to_mode_xn,to_mode_num,status,dxp] = transition(obj,from_mode_num,t,mode_x,u)
                mode_x = double(mode_x);
                u = double(u);
                dt = 0.05;
                [xdot,dxdot] = impact_dynamics(obj,t,mode_x,u,dt);
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
        
        function [xdot, dxdot] = impact_dynamics(obj,t,x,u,dt)
            F_impact = 0;
            if nargout == 1
                %xdot = impact_foamy_pendulum_dynamics(t,x,u,obj.m_package,F_impact,dt);
                xdot = impact_foamy_pendulum_dynamics_mex(t,x,u,obj.m_package,dt);
                dxdot = 0;
            else %Need Jacobian
                options.grad_level = 0;
                options.grad_method = 'numerical';
                options.diff_type = 'central';
                %[xdot, dxdot] = geval(@(t1,x1,u1) impact_foamy_pendulum_dynamics(t1,x1,u1,obj.m_package,F_impact,dt),t,x,u,options);
                [xdot, dxdot] = geval(@(t1,x1,u1) impact_foamy_pendulum_dynamics_mex(t1,x1,u1,obj.m_package,dt),t,x,u,options);
            end
        end
        
        function y = output(obj,t,x,u)
            y = x;
        end
        
        function [xtraj,utraj] = runDircol(obj,display)

            % initial conditions:
            [x0, u0] = findTrim(obj,5); %find trim conditions for level flight at 6 m/s
            x0(1) = -6;
            x0(3) = x0(3) + 1.5;
            x0(8) = x0(8) - 6;
            x0(10) = x0(10)+1.5;
            
            % final conditions:
            xf = x0;
            xf(1) = x0(1)+12; %translated in x
            xf(8) = x0(8)+12;

            tf0 = (xf(1)-x0(1))/6; % initial guess at duration 

            N = [17;17];
            
            options.u_const_across_transitions = false;
            
            prog = HybridTrajectoryOptimization(@DircolTrajectoryOptimization,obj,[1;2],N,{[0,tf0/2];[0,tf0/2]},options);
            %prog = DircolTrajectoryOptimization(obj,N,[0 tf0]);
            prog = addModeStateConstraint(prog,1,ConstantConstraint(x0),1);
            prog = addModeStateConstraint(prog,2,ConstantConstraint(xf),N(2));
            prog = addModeStateConstraint(prog,1,QuadraticConstraint(.5,.5,eye(4),zeros(4,1)),1:N(1),4:7);
            prog = addModeStateConstraint(prog,2,QuadraticConstraint(.5,.5,eye(4),zeros(4,1)),1:N(2),4:7);
            prog = addModeInputConstraint(prog,1,BoundingBoxConstraint([0.001; -1; -1; -1], [1; 1; 1; 1]),1:N(1));
            prog = addModeInputConstraint(prog,2,BoundingBoxConstraint([0.001; -1; -1; -1], [1; 1; 1; 1]),1:N(2));
            prog = addModeRunningCost(prog,1,@cost);
            prog = addModeRunningCost(prog,2,@cost);
            %prog = addFinalCost(prog,@(t,x) finalCost(t,x,xf));

            %--- snopt options ---%

            %--- snopt options ---%
            prog = setSolver(prog,'snopt');
            prog = prog.setSolverOptions('snopt','minoroptimalitytolerance',1e-4);
            prog = prog.setSolverOptions('snopt','minorfeasibilitytolerance',1e-4);
            prog = prog.setSolverOptions('snopt','majoroptimalitytolerance',1e-4);
            prog = prog.setSolverOptions('snopt','majorfeasibilitytolerance',1e-4);
            prog = prog.setSolverOptions('snopt','iterationslimit',100000);

            t_init_full = linspace(0,tf0,N(1)+N(2));
            t_init{1} = t_init_full(1:N(1));
            t_init{2} = t_init_full(N(1)+1:end);
            
            %Warm Starting
            %[xtraj_temp,utraj_temp] = obj.modes{1}.runDircol(0);

            %Set initial guess for controls to be trim conditions
            traj_init{1}.u = setOutputFrame(PPTrajectory(foh(t_init{1},kron(ones(1,N(1)),u0))),getInputFrame(obj));
            traj_init{2}.u = setOutputFrame(PPTrajectory(foh(t_init{2},kron(ones(1,N(2)),u0))),getInputFrame(obj));
         
            
            %Simulate with u0 input to generate initial guess
            [t_guess1, x_guess1] = ode45(@(t,x) obj.dynamics(t,[1;x],u0),t_init{1},x0);
            [t_guess2, x_guess2] = ode45(@(t,x) obj.dynamics(t,[2;x],u0),t_init{2},x_guess1(end,:));
            t_guess{1} = t_guess1;
            t_guess{2} = t_guess2;
            x_guess{1} = x_guess1;
            x_guess{2} = x_guess2;
            %t_guess = [t_guess;t_guess2];
            %x_guess = [x_guess;x_guess2];
            %x_guess = [[ones(5,1);2*ones(5,1)],x_guess]
            
            traj_init{1}.x = setOutputFrame(PPTrajectory(foh(t_guess{1},x_guess{1}')),getOutputFrame(obj));
            traj_init{2}.x = setOutputFrame(PPTrajectory(foh(t_guess{2},x_guess{2}')),getOutputFrame(obj));
            

            
            
            %traj_init{1}.x0 = x0';
            %traj_init{2}.x0 = x_guess1(end,:)';
            
            prog = prog.compile();
            tic
            [xtraj,utraj,~,~,info]=solveTraj(prog,t_init,traj_init);
            %[xtraj,utraj,~,~,info]=solveTraj(prog,[t_init/2;t_init/2]);
            toc


            if nargin == 2 && display
               visualizeFoamy(obj,xtraj,true);
            end

            function [g,dg] = cost(dt,x,u)
                R = eye(4);
                g = 0.5*(u-u0)'*R*(u-u0);
                dg = [0, zeros(1,26), (u-u0)'*R];
            end

            function [h,dh] = finalCost(t,x,xf)
                hx = .5*(x(1:3)-xf(1:3))'*(x(1:3)-xf(1:3));
                %hv = .5*(x(8:10)-xf(8:10))'*(x(8:10)-xf(8:10));
                %hw = .5*(x(11:13)-xf(11:13))'*(x(11:13)-xf(11:13));

                %This is to handle the quaternion double cover issue
                hq1 = 1 - xf(4:7)'*x(4:7);
                hq2 = 1 + xf(4:7)'*x(4:7);
                if hq1 <= hq2
                    h = hx + hq1;
                    dh = [0,x(1:3)'-xf(1:3)',-xf(4:7)',zeros(1,6)];
                else
                    h = hx + hq2;
                    dh = [0,x(1:3)'-xf(1:3)',xf(4:7)',zeros(1,6)];
                end
            end
            
        end
        
        function [xtrim, utrim] = findTrim(obj,v,varargin)
            %Trim the plane for level forward flight at speed v

            v0 = [v 0 0]';
            w0 = [0 0 0]';
            q0 = [0 1 0 0]';
            u0 = [.5; 0; 0; 0; 0; 0; 0; 0; -.35];

            function r = trimResidual(x)
                if isempty(varargin)
                    q = [x(5) 1 x(6) x(7)]';
                else
                    q = [1 x(5) x(6) x(7)]';
                end
                q = q/sqrt(q'*q);
                y = [1;0;0;0;q;x(8);0;x(9);q0;v0;w0;v0;w0];
                xdot = obj.dynamics(0,y,x(1:4));
                r = [xdot(4:7); xdot(11:26)]; 
            end

            options = optimoptions('fsolve','Algorithm','levenberg-marquardt','FiniteDifferenceType','central','Display','off','TolFun',1e-10,'TolX',1e-10);
            y = fsolve(@trimResidual,u0,options);

            if isempty(varargin)
                qtrim = [y(5) 1 y(6) y(7)]';
            else
                qtrim = [1 y(5) y(6) y(7)]';
            end
            qtrim = qtrim/sqrt(qtrim'*qtrim);
            xtrim = [0;0;0;qtrim;y(8);0;y(9);q0;v0;w0;v0;w0];
            utrim = y(1:4);

        end

        function x = getInitialState(obj)
            x = zeros(26,1);
            x(5) = 1; %Quaternion - plane oriented right-side up
            x(10) = -.35; %Pendulum position under airplane
            x(12) = 1; %Pendulum quaternion
        end
        
        function q = getZeroConfiguration(obj)
            q = zeros(14,1);
            q(5) = 1; %Quaternion - plane oriented right-side up
            q(10) = -.35; %Pendulum position under airplane
            q(12) = 1; %Pendulum quaternion
        end
    end
    
end



    

