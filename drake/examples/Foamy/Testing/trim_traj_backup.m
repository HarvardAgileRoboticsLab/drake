       function [xtraj,utraj] = runDircol_trim(obj,display)

            % initial conditions:
            [x0, u0] = findTrim(obj,5); %find trim conditions for level flight at 6 m/s
            x0(1) = -10;
            x0(3) = x0(3) + 1.5;
            x0(8) = x0(8) - 10;
            x0(10) = x0(10)+1.5;
            
            % final conditions:
            xf = x0;
            xf(1) = x0(1)+20; %translated in x
            xf(8) = x0(8)+20;

            tf0 = (xf(1)-x0(1))/6; % initial guess at duration 

            N = [15;15];
            
            N_trim_start = 12; %Index to start trim condition
            trim_length = 3; %Length, in knot points, of trim condition
            trim_inds = {};
        
            for i = 1:trim_length
                trim_inds{i} = [N_trim_start-1+i,N_trim_start+i];   
            end

            
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
            
            %num_trims = [4 5 6 7 15 16 17 18 19 20];
            trim_vars = [15 16 17 18 19 20];
            num_trims = length(trim_vars);
%             for j = 1:length(num_trims)
%                 prog = addModeStateConstraint(prog,1,FunctionHandleConstraint(0,0,2,@(x1,x2) trim_constraint(x1,x2,j)),trim_inds,num_trims(j));
%             end
            
            prog = addModeStateConstraint(prog,1,FunctionHandleConstraint(zeros(num_trims,1),zeros(num_trims,1),num_trims*2,@trim_constraint),trim_inds,trim_vars);
            
            %prog = addFinalCost(prog,@(t,x) finalCost(t,x,xf));


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
            
         function [c,dc] = trim_constraint(x1,x2)
             l = length(x1);
             c = (x2-x1);
             dc = [-1*eye(l),eye(l)];
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