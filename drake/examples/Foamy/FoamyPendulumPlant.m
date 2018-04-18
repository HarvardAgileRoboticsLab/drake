classdef FoamyPendulumPlant < DrakeSystem
    
    properties
        m_package = 0; %Package mass
        disturbance_type = 1; %1 = package mass, 2 = airplane mass, 3 = package dynamcis, 4 = wind gusts
       
        
    end
    
    methods
        
        function obj = FoamyPendulumPlant(m_p)
            obj = obj@DrakeSystem(26,0,4,26,0,1);
            obj = obj.setStateFrame(CoordinateFrame('FoamyPendulumState',26,'x',{'x1','x2','x3','q0','q1','q2','q3','x4','x5','x6','q4','q5','q6','q7','v1','v2','v3','w1','w2','w3','v4','v5','v6','w4','w5','w6'}));
            obj = obj.setInputFrame(CoordinateFrame('FoamyPendulumInput',4,'u',{'thr','ail','elev','rud'}));
            obj = obj.setOutputFrame(obj.getStateFrame); %full state feedback
            %quaternion unit-norm constraint -- why doesn't the index input work?
            %obj = obj.addStateConstraint(QuadraticConstraint(.5,.5,blkdiag(zeros(3),eye(4),zeros(6)),zeros(13,1)));
            if nargin == 0
                obj.m_package = 0;
            else
                obj.m_package = m_p;
            end
        end
           
           
   function nw = getNumDisturbances(obj)
        switch obj.disturbance_type
            case 1
                nw = 1;
            case 2
                nw = 1;
            case 3
                nw = obj.getNumContStates();
            case 4
                nw = 3;
            otherwise
                error('Unknown disturbance type');
        end
    end

        
        function obj = change_m_package(obj,m_p)
            obj = set.m_package(obj,m_p);
       % obj = obj.m_package = m_p;
        end

        function [xdot, dxdot] = dynamics(obj,t,x,u)
            if nargout == 1
                %xdot = foamy_pendulum_dynamics(t,x,u);
                xdot = foamy_pendulum_dynamics_mex(t,x,u,obj.m_package);
                dxdot = 0;
            else %Need Jacobian
                options.grad_level = 0;
                options.grad_method = 'numerical';
                options.diff_type = 'central';
                %[xdot, dxdot] = geval(@(t1,x1,u1) foamy_pendulum_dynamics(t1,x1,u1),t,x,u,options);
                [xdot, dxdot] = geval(@(t1,x1,u1) foamy_pendulum_dynamics_mex(t1,x1,u1,obj.m_package),t,x,u,options);
            end
        end
        
  function [f,df,d2f] = dynamics_w(obj,t,x,u,w)
        
        [f,df] = dynamics_w_mp(obj,t,x,u,w);
        
        if nargout == 3
            %Finite diff to get 2nd derivatives
            nx = length(x);
            nu = length(u);
            nw = length(w);
            
            Dx = 1e-6*eye(nx);
            Du = 1e-6*eye(nu);
            Dw = 1e-6*eye(nw);
            
            d2f = zeros(nx, 1+nx+nu+nw, 1+nx+nu+nw);
            for k = 1:nx
                [~,df_p] = dynamics_w_mp(obj,t,x+Dx(:,k),u,w);
                d2f(:,:,1+k) = df_p-df;
            end
            for k = 1:nu
                [~,df_p] = dynamics_w_mp(obj,t,x,u+Du(:,k),w);
                d2f(:,:,1+nx+k) = df_p-df;
            end
            for k = 1:nw
                [~,df_p] = dynamics_w_mp(obj,t,x,u,w+Dw(:,k));
                d2f(:,:,1+nx+nu+k) = df_p-df;
            end
            
            d2f = reshape(d2f,nx,(1+nx+nu+nw)*(1+nx+nu+nw));
            
        end

  end
    
 function [f,df] = dynamics_w_mp(obj,t,x,u,w)
      % w is the weight difference for the perceived package

      nq = obj.getNumPositions;
      q=x(1:nq); 
      
      if (nargout>1)
        nx = obj.getNumStates;
        nu = obj.getNumInputs;

        %kinsol = doKinematics(obj, q, [], struct('compute_gradients', true));
        %[~,J,dJ] = obj.forwardKin(kinsol,findLinkId(obj,obj.hand_name),[0;0;0]);
        %uw = J'*w;
 
        %dJtw = zeros(nq,nq);
        %for i=1:nq
        %  dJtw(:,i) = dJ(:,(i-1)*nq+(1:nq))'*w;
        %end
      
        [f,df] = geval(@(t1,x1,u1,w1) foamy_pendulum_dynamics_mex(t1,x1,u1,obj.m_package+w1),t,x,u,w,options);
        df_du = df(:,1+nx+(1:nu)); 
        df_dq = df(:,1+(1:nq));
        %df_dq = df(:,1+(1:nq)) + df_du*dJtw;
        df_dqd = df(:,1+nq+(1:nq));
        %df_dw = df_du*J';
        df_dw = df(:,1+nq+(1:nq));
        
        df = [df(:,1),df_dq,df_dqd,df_du,df_dw];
      else
        %kinsol = doKinematics(obj, q, []);
        %[~,J] = obj.forwardKin(kinsol,findLinkId(obj,obj.hand_name),[0;0;0]);
        %uw = J'*w;
      
        %[f,df] = dynamics(obj,t,x,u+uw);
        
        [f,df] = geval(@(t1,x1,u1,w1) foamy_pendulum_dynamics_mex(t1,x1,u1,obj.m_package+w1),t,x,u,w,options);
      end

    end
        
        function y = output(obj,t,x,u)
            y = x;
        end
        
        function [xtraj,utraj] = runDircol(obj,display)

            % initial conditions:
            [x0, u0] = findTrim(obj,6); %find trim conditions for level flight at 8 m/s
            x0(1) = x0(1)-6;
            x0(3) = x0(3)+1.5;
            x0(8) = x0(8)-6;
            x0(10) = x0(10)+1.5;

            % final conditions:
            xf = x0;
            xf(1) = x0(1)+12; %translated in x
            xf(8) = x0(8)+12;

            tf0 = (xf(1)-x0(1))/6; % initial guess at duration 

            N = 20;
            prog = DircolTrajectoryOptimization(obj,N,[0 tf0]);
            prog = addStateConstraint(prog,ConstantConstraint(x0),1);
            prog = addStateConstraint(prog,ConstantConstraint(xf),N);
            prog = addStateConstraint(prog,QuadraticConstraint(.5,.5,eye(4),zeros(4,1)),1:N,4:7);
            prog = addStateConstraint(prog,QuadraticConstraint(.5,.5,eye(4),zeros(4,1)),1:N,11:14);
            prog = addInputConstraint(prog,BoundingBoxConstraint([.001; -1; -1; -1], [1; 1; 1; 1]),1:N);
            prog = addRunningCost(prog,@cost);
            prog = addStateConstraint(prog,ConstantConstraint(0.5),10,10);
            %prog = addFinalCost(prog,@(t,x) finalCost(t,x,xf));

            %--- snopt options ---%
            prog = setSolver(prog,'snopt');
            prog = prog.setSolverOptions('snopt','minoroptimalitytolerance',1e-4);
            prog = prog.setSolverOptions('snopt','minorfeasibilitytolerance',1e-4);
            prog = prog.setSolverOptions('snopt','majoroptimalitytolerance',1e-4);
            prog = prog.setSolverOptions('snopt','majorfeasibilitytolerance',1e-4);
            prog = prog.setSolverOptions('snopt','iterationslimit',100000);

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
                v = FoamyPendulumVisualizer(obj);
                v.playback(xtraj,struct('slider',true));
            end

            function [g,dg] = cost(dt,x,u)
                R = eye(4);
                g = 0.5*(u-u0)'*R*(u-u0);
                dg = [0, zeros(1,26), (u-u0)'*R];
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
                y = [0;0;0;q;x(8);0;x(9);q0;v0;w0;v0;w0];
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

