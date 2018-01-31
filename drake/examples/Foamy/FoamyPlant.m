classdef FoamyPlant < DrakeSystem
    
    properties
        
    end
    
    methods
        
        function obj = FoamyPlant()
            obj = obj@DrakeSystem(13,0,4,13,0,1);
            obj = obj.setStateFrame(CoordinateFrame('FoamyState',13,'x',{'x1','x2','x3','q0','q1','q2','q3','v1','v2','v3','w1','w2','w3'}));
            obj = obj.setInputFrame(CoordinateFrame('FoamyInput',4,'u',{'thr','ail','elev','rud'}));
            obj = obj.setOutputFrame(obj.getStateFrame); %full state feedback
            %quaternion unit-norm constraint -- why doesn't the index input work?
            %obj = obj.addStateConstraint(QuadraticConstraint(.5,.5,blkdiag(zeros(3),eye(4),zeros(6)),zeros(13,1)));
        end
        
        function [xdot, dxdot] = dynamics(obj,t,x,u)
            if nargout == 1
                %xdot = foamy_dynamics(t,x,u);
                xdot = foamy_dynamics_mex(t,x,u);
                dxdot = 0;
            else %Need Jacobian
                options.grad_level = 0;
                options.grad_method = 'numerical';
                options.diff_type = 'central';
                %[xdot, dxdot] = geval(@(t1,x1,u1) foamy_dynamics(t1,x1,u1),t,x,u,options);
                [xdot, dxdot] = geval(@(t1,x1,u1) foamy_dynamics_mex(t1,x1,u1),t,x,u,options);
            end
        end
        
        function [xdot, dxdot,d2xdot] = dynamics_w(obj,t,x,u,w)
            if nargout == 1
                xdot = varMass_foamy_dynamics(t,x,u,w);
                %xdot = foamy_dynamics_mex(t,x,u);
                dxdot = 0;
                d2xdot = 0;
            else %Need Jacobian
                options.grad_level = 0;
                options.grad_method = 'numerical';
                options.diff_type = 'central';
                [xdot, dxdot] = geval(@(t1,x1,u1) varMass_foamy_dynamics(t1,x1,u1,w),t,x,u,options);
                d2xdot = 0;
                %[xdot, dxdot] = geval(@(t1,x1,u1) foamy_dynamics_mex(t1,x1,u1),t,x,u,options);
            end
        end
        
        function y = output(obj,t,x,u)
            y = x;
        end
        
        function n = getNumDisturbances(~)
            n=1;
        end
        
        function [xtraj,utraj] = runDircol(obj,display)

            % initial conditions:
            [x0, u0] = findTrim(obj,6); %find trim conditions for level flight at 6 m/s
            x0(1) = -6;
            x0(3) = 1.5;

            % final conditions:
            xf = x0;
            xf(1) = 6; %translated in x

            tf0 = (xf(1)-x0(1))/6; % initial guess at duration 

            N = 7;
            prog = DircolTrajectoryOptimization(obj,N,[0 tf0]);
            prog = addStateConstraint(prog,ConstantConstraint(x0),1);
            prog = addStateConstraint(prog,ConstantConstraint(xf),N);
            prog = addStateConstraint(prog,QuadraticConstraint(.5,.5,eye(4),zeros(4,1)),1:N,4:7);
            prog = addInputConstraint(prog,BoundingBoxConstraint([0; -1; -1; -1], [1; 1; 1; 1]),1:N);
            prog = addRunningCost(prog,@cost);
            %prog = addFinalCost(prog,@(t,x) finalCost(t,x,xf));

            %--- snopt options ---%
            prog = setSolver(prog,'snopt');
            prog = prog.setSolverOptions('snopt','majoroptimalitytolerance',1e-5);
            prog = prog.setSolverOptions('snopt','majorfeasibilitytolerance',1e-5);

            t_init = linspace(0,tf0,N);

            %Set initial guess for controls to be trim conditions
            traj_init.u = setOutputFrame(PPTrajectory(foh(t_init,kron(ones(1,N),u0))),getInputFrame(obj));

            %Simulate with u0 input to generate initial guess
            [t_guess, x_guess] = ode45(@(t,x) obj.dynamics(t,x,u0),t_init,x0);
            traj_init.x = setOutputFrame(PPTrajectory(foh(t_guess,x_guess')),getStateFrame(obj));

            tic
            [xtraj,utraj,~,~,info]=solveTraj(prog,t_init,traj_init);
            toc

            if nargin == 2 && display
               visualizeFoamy(obj,xtraj,true);
            end

            function [g,dg] = cost(dt,x,u)
                R = eye(4);
                g = 0.5*(u-u0)'*R*(u-u0);
                dg = [0, zeros(1,13), (u-u0)'*R];
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
        
     function [utraj,xtraj,z,prog] = robustTrajectory(obj,N,D,options)
        
        if nargin == 3
            options = struct();
        end
        
        % initial conditions:
        [x0, u0] = findTrim(obj,6); %find trim conditions for level flight at 6 m/s
        x0(1) = -6;
        x0(3) = 1.5;

        % final conditions:
        xf = x0;
        xf(1) = 6; %translated in x
        
        tf0 = (xf(1)-x0(1))/6; % initial guess at duration 
        
        Q = eye(13);
        Q(1,1) = 10;
        R = eye(4);
        Qf = 100*eye(13);
        
        E0 = zeros(13,13);
        
        prog = RobustDirtranTrajectoryOptimization(obj,N,D,E0,Q,R,Qf,[0 tf0],options);
        prog = prog.addStateConstraint(ConstantConstraint(x0),1);
        prog = prog.addStateConstraint(ConstantConstraint(xf),N);
        prog = addStateConstraint(prog,QuadraticConstraint(.5,.5,eye(4),zeros(4,1)),1:N,4:7);
        %prog = addInputConstraint(prog,BoundingBoxConstraint([0; -1; -1; -1], [1; 1; 1; 1]),1:N);
        prog = addRunningCost(prog,@cost);
        %prog = prog.addFinalCost(@finalCost);
        
        prog = prog.addRobustCost(Q,R,Qf);
        prog = prog.addRobustInputConstraint();
        
        prog = prog.setSolverOptions('snopt','majoroptimalitytolerance', 1e-2);
        prog = prog.setSolverOptions('snopt','majorfeaasibilitytolerance', 1e-3);
        prog = prog.setSolverOptions('snopt','minorfeaasibilitytolerance', 1e-3);
        
        %prog = addTrajectoryDisplayFunction(prog,@displayTrajectory);
        
        traj_init.x = PPTrajectory(foh([0,tf0],[double(x0),double(xf)]));
        
        disp('Running solve');
        tic
        [xtraj,utraj,z] = prog.solveTraj(tf0,traj_init);
        toc
        
        function [h,dh] = finalCost(tf,x)
            h = tf;
            if (nargout>1)
                dh = [1, zeros(1,13)];
            end
        end
        function [g,dg] = cost(dt,x,u)
                R = eye(4);
                g = 0.5*(u-u0)'*R*(u-u0);
                dg = [0, zeros(1,13), (u-u0)'*R];
        end
        
    end
        
        function [xtrim, utrim] = findTrim(obj,v,varargin)
            %Trim the plane for level forward flight at speed v

            v0 = [v 0 0]';
            w0 = [0 0 0]';

            u0 = [.5; 0; 0; 0; 0; 0; 0];

            function r = trimResidual(x)
                if isempty(varargin)
                    q = [x(5) 1 x(6) x(7)]';
                else
                    q = [1 x(5) x(6) x(7)]';
                end
                q = q/sqrt(q'*q);
                y = [0;0;0;q;v0;w0];
                xdot = obj.dynamics(0,y,x(1:4));
                r = xdot(4:13);
            end

            options = optimoptions('fsolve','Algorithm','levenberg-marquardt','FiniteDifferenceType','central','Display','off','TolFun',1e-10,'TolX',1e-10);
            y = fsolve(@trimResidual,u0,options);

            if isempty(varargin)
                qtrim = [y(5) 1 y(6) y(7)]';
            else
                qtrim = [1 y(5) y(6) y(7)]';
            end
            qtrim = qtrim/sqrt(qtrim'*qtrim);
            xtrim = [0;0;0;qtrim;v0;w0];
            utrim = y(1:4);

        end

        function x = getInitialState(obj)
            x = zeros(13,1);
            x(5) = 1; %Quaternion - plane oriented right-side up
        end
        
        function q = getZeroConfiguration(obj)
            q = zeros(7,1);
            q(5) = 1; %Quaternion - plane oriented right-side up
        end
    end
    
end

