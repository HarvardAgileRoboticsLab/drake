classdef McfoamyModel < DrakeSystem
    
    properties
        
    end
    
    methods
        
        function obj = McfoamyModel()
            obj = obj@DrakeSystem(13,0,4,13,0,1);
            obj = obj.setStateFrame(CoordinateFrame('ExtraState',13,'x',{'x1','x2','x3','q0','q1','q2','q3','v1','v2','v3','w1','w2','w3'}));
            obj = obj.setInputFrame(CoordinateFrame('ExtraInput',4,'u',{'thr','ail','elev','rud'}));
            obj = obj.setOutputFrame(obj.getStateFrame); %full state feedback
            %quaternion unit-norm constraint -- why doesn't the index input work?
            obj = obj.addStateConstraint(QuadraticConstraint(.5,.5,blkdiag(zeros(3),eye(4),zeros(6)),zeros(13,1)));
        end
        
        function [xdot, dxdot] = dynamics(obj,t,x,u)
            if nargout == 1
                xdot = foamy_plant(t,x,u);
                %xdot = extra_dynamics_mex(t,x,u);
                dxdot = 0;
            else %Need Jacobian
                options.grad_level = 0;
                options.grad_method = 'numerical';
                options.diff_type = 'central';
                [xdot, dxdot] = geval(@(t1,x1,u1) foamy_plant(t1,x1,u1),t,x,u,options);
                %[xdot, dxdot] = geval(@(t1,x1,u1) extra_dynamics_mex(t1,x1,u1),t,x,u,options);
            end
        end
        
        function y = output(obj,t,x,u)
            y = x;
        end
        
        function controller = prophang(obj, display)
            %Simulate an LQR-stabilized prophang
            x0 = [0 0 1.5 1/sqrt(2) 0 -1/sqrt(2) 0 0 0 0 0 0 0]';
            u0 = [135.893015369865 106 106 106]';
            Q = blkdiag(eye(3), .1*eye(4), 10*eye(6));
            R = 1e-6*eye(4);
            
            ltisys = tilqr(obj,x0,u0,Q,R);
            clsys = feedback(obj,ltisys);
            
            %Not sure why this is necessary, but the numbers that come out
            %of the controller are weird if it's not first made into a
            %feedback system with the plant.
            controller = clsys.sys2;
            
            %Plot
            if nargin == 2 && display
                %Perturb initial position
                xsim = x0;
                xsim(1:3) = xsim(1:3)+[.4; .1; -.3];

                %Simulate closed loop system
                [~,cltraj] = clsys.simulateODE([0 10],xsim);
            
                v = ExtraVisualizer(obj);
                v.playback(cltraj);
            end
        end
        
        function lqgsystem = prophang2(obj, display)
            %Simulate an LQR-stabilized prophang with observer
            [x0, u0] = findHover(obj);
            x0(1:3) = [.5 0 1]';
            
            q0 = [.72; -.05; -.69; -.04]; %Missalignment of body frame axes
            q0 = q0/norm(q0);
            
            [A,B] = obj.linearize(0,x0,u0);
            P = blkdiag(eye(3),[-q0(2:4), q0(1)*eye(3) - hat(q0(2:4))], eye(6));
            A = P*A*P';
            B = P*B;
            
            u0 = [147; 39; 106; 100]; %Poorly tuned u0
            
            Q = blkdiag(diag([5 5 10]), diag([300 2 2]), 30*eye(3), 30*eye(3));
            R = diag([1e-4, 1e-4, 1e-4, 1e-4]);
            K = lqr(A,B,Q,R);
            %controller = AffineSystem([],[],[],[],[],[],[],-K*P,u0+K*P*x0);
            %controller = controller.setInputFrame(obj.getOutputFrame);
            %controller = controller.setOutputFrame(obj.getInputFrame);
            %clsys = feedback(obj,controller);
            
            trimmer = LinearSystem([],[],[],[],[],[eye(7) zeros(7,6)]);
            trimmer = trimmer.setInputFrame(obj.getOutputFrame);
            plant = cascade(obj,trimmer);
            
            C = [eye(6) zeros(6)];
            sys = ss(A,[B eye(12)],C,[zeros(6,4) zeros(6,12)]);
            Qe = eye(12);
            Re = 1e-4*eye(6);
            [~,L] = kalman(sys,Qe,Re);
            
            estimator = AffineSystem((A-B*K-L*C),L*P(1:6,1:7),-L*P(1:6,1:7)*x0(1:7),[],[],[],eye(12),zeros(12,7),zeros(12,1));
            controller = AffineSystem([],[],[],[],[],[],[],-K,u0);
            controller = controller.setInputFrame(estimator.getOutputFrame);
            lqgsystem = cascade(estimator,controller);
            lqgsystem = lqgsystem.setInputFrame(plant.getOutputFrame);
            lqgsystem = lqgsystem.setOutputFrame(plant.getInputFrame);
            
            clsys = feedback(plant,lqgsystem);
            
            %Plot
            if nargin == 2 && display
                %Perturb initial position
                xsim = x0;
                xsim(1:3) = xsim(1:3)+[.04; .02; -.03];
                xsim = [xsim; zeros(12,1)];

                %Simulate closed loop system
                [~,cltraj] = clsys.simulateODE([0 10],xsim);
                extratraj = cltraj(1:13);
                extratraj = extratraj.setOutputFrame(obj.getStateFrame);
                v = ExtraVisualizer(obj);
                v.playback(extratraj);
            end
        end
        
        function [cltraj, controller] = runLQR(obj,display)

            %Generate nominal trajectory
            [xtraj,utraj] = runDircol(obj);

            %Weighting matrices for TVLQR
            Q = diag([(1/0.05)^2 (1/0.05)^2 (1/0.05)^2 .1 .1 .1 .1 1 1 1 1 1 1]);
            Qf = diag([(1/0.05)^2 (1/0.05)^2 (1/0.05)^2 (1/0.05)^2 (1/0.05)^2 (1/0.05)^2 (1/0.05)^2 1 1 1 1 1 1]);
            R = diag([.1 .1 .1 .1]);

            ltvsys = tvlqr(obj,xtraj,utraj,Q,R,Qf);

            clsys = feedback(obj,ltvsys);

            %Not sure why this is necessary, but the numbers that come out
            %of the controller are weird if it's not first made into a
            %feedback system with the plant.
            controller = clsys.sys2;

            %Perturb initial position and velocity
            x0 = xtraj.eval(0);
            x0(1:3) = [-2.6; -.5;.5];
            x0(8:10) = [5.5; -.1; .1];

            %Simulate closed loop system
            [~,cltraj] = clsys.simulateODE([0 1],x0);

            %Plot
            if nargin == 2 && display
                v = ExtraVisualizer(obj);
                v.playback(cltraj);
            end

        end
        
        function [xtraj,utraj]=runDircol(obj,display)

            % initial conditions:
            [x0, u0] = findTrim(obj,7); %find trim conditions for level flight at 6 m/s
            x0(1) = -2.5;
            x0(3) = 1.5;

            % final conditions:
            xf = x0;
            xf(1) = 3.5; %translated in x

            tf0 = 6/7; % initial guess at duration 

            N = 7;
            prog = DircolTrajectoryOptimization(obj,N,[0 tf0]);
            prog = addStateConstraint(prog,BoundingBoxConstraint(x0-1e-4,x0+1e-4),1); %Relax the initial condition constraint a little
            prog = addRunningCost(prog,@cost);
            prog = addFinalCost(prog,@(t,x) finalCost(t,x,xf));

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
                v = ExtraVisualizer(obj);
                v.playback(xtraj);
            end

            function [g,dg] = cost(dt,x,u)
                g = 0;
                dg = zeros(1, 1 + size(x,1) + size(u,1));
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

            u0 = [40; 90; .04];

            function c = cost(x)
                if isempty(varargin)
                    q = [0 1 0 x(3)]';
                else
                    q = [1 0 x(3) 0]';
                end
                q = q/sqrt(q'*q);
                y = [0;0;0;q;v0;w0];
                xdot = obj.dynamics(0,y,[x(1); 107; x(2); 107]);
                c = sqrt(xdot(4:13)'*xdot(4:13));
            end

            options = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton','TolFun',1e-10,'TolX',1e-10);
            y = fminunc(@cost,u0,options);

            if isempty(varargin)
                qtrim = [0 1 0 y(3)]';
            else
                qtrim = [1 0 y(3) 0]';
            end
            qtrim = qtrim/sqrt(qtrim'*qtrim);
            xtrim = [0;0;0;qtrim;v0;w0];
            utrim = [y(1); 106; y(2); 106];

        end
        
        function [xhover, uhover] = findHover(obj)
            %find control inputs for hover
            x0 = [0 0 0 1/sqrt(2) 0 -1/sqrt(2) 0 0 0 0 0 0 0]';
            u0 = [135 35 107 107]';
            
            function c = cost(x)
                xdot = obj.dynamics(0,x0,[x; u0(3:4)]);
                c = sqrt(xdot'*xdot);
            end
            
            options = optimoptions(@fminunc,'Display','off','Algorithm','quasi-newton','TolFun',1e-10,'TolX',1e-10);
            y = fminunc(@cost,u0(1:2),options);
            
            xhover = x0;
            uhover = [y; u0(3:4)];
        end

        function x = getInitialState(obj)
            x = zeros(13,1);
            x(5) = 1; %Quaternion - plane oriented right-side up
        end
        
    end
    
end



