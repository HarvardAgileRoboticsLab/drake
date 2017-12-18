classdef CUDPTrajectoryOptimization < IndirectTrajectoryOptimization
    %Iterative LQR Trajectory Optimization + Constraints
    
    properties (Constant)
        EULER = 1;
        MIDPOINT = 2;
        RK3 = 3; % DEFAULT
        RK4 = 4;
        iLQR = 1;
        UDP = 2;
    end
    
    properties (Access = private)
        integrator
        integrator_b
        numXinConstraints
        numXeqConstraints
        numUinConstraints
        numUeqConstraints
        xinConstraint
        xeqConstraint
        uinConstraint
        ueqConstraint
        planVisualizer
        computeReceedingMultiple
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTRUCTOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = CUDPTrajectoryOptimization(plant,N,duration,options)
            
            defaults = struct(...
                'maxMinorIter',     1000,...   maximum iterations
                'maxMajorIter',     100,...     maximum iterations
                'rhoFactor',        4,...       rho scaling factor
                'rhoMax',           1e5,...     maximum rho value
                'rhoMin',           1e-5,...    minimum rho value
                'alphaMin',         1e-4,...    minimum line search parameter
                'alphaFactor',      2,...       alpha scaling factor
                'muInit',           1,...       initial mu value
                'muFactor',         10,...      mu scaling factor
                'muMax',            1e12,...    maximum mu value
                'tolFun',           1e-2,...    minimum cost reduction
                'tolFunExit',       1e-6,...    minimum cost reduction on final iter
                'tolGrad',          1e-1,...    minimum gradient size
                'tolGradExit',      1e-4,...    minimum gradient size on final iter
                'tolConstr',        1e-5,...    constraint tolerance for exit
                'tolLambdaUpdate',  5e-2,...    constraint violation minimum for lambda updates
                'phiFactor',        10,...      factor for phi (lambda update tolerance)
                'ILambdaDecay',     1,...       number of major iterations for decay of the ILambda term for inequality constraints
                'penaltyOnly',      false,...   only use the penalty term no lambda updates
                'BoxQPActive',      false,...   whether or not to run the BoxQP (only works with iLQR)
                'BoxQPBounds',      [],...      [Lower,Upper] bounds for inputs (only works with iLQR)
                'SigmaScale',       1,...       Scaling factor for sigma points
                'receedingHorizon', false,...   receeding horizon on tolerance
                'debug',            false,...   print debug data?
                'derivative_method', CUDPTrajectoryOptimization.UDP,... iLQR vs UDP
                'integration_method', CUDPTrajectoryOptimization.RK3);
            
            if nargin < 4
                options = defaults;
            end
            if ~isfield(options,'maxMinorIter')
                options.maxMinorIter = defaults.maxMinorIter;
            end
            if ~isfield(options,'maxMajorIter')
                options.maxMajorIter = defaults.maxMajorIter;
            end
            if ~isfield(options,'rhoFactor')
                options.rhoFactor = defaults.rhoFactor;
            end
            if ~isfield(options,'rhoMax')
                options.rhoMax = defaults.rhoMax;
            end
            if ~isfield(options,'rhoMin')
                options.rhoMin = defaults.rhoMin;
            end
            if ~isfield(options,'alphaMin')
                options.alphaMin = defaults.alphaMin;
            end
            if ~isfield(options,'alphaFactor')
                options.alphaFactor = defaults.alphaFactor;
            end
            if ~isfield(options,'muInit')
                options.muInit = defaults.muInit;
            end
            if ~isfield(options,'muFactor')
                options.muFactor = defaults.muFactor;
            end
            if ~isfield(options,'muMax')
                options.muMax = defaults.muMax;
            end
            if ~isfield(options,'tolFun')
                options.tolFun = defaults.tolFun;
            end
            if ~isfield(options,'tolFunExit')
                options.tolFunExit = defaults.tolFunExit;
            end
            if ~isfield(options,'tolGrad')
                options.tolGrad = defaults.tolGrad;
            end
            if ~isfield(options,'tolGradExit')
                options.tolGradExit = defaults.tolGradExit;
            end
            if ~isfield(options,'tolConstr')
                options.tolConstr = defaults.tolConstr;
            end
            if ~isfield(options,'tolLambdaUpdate')
                options.tolLambdaUpdate = defaults.tolLambdaUpdate;
            end
            if ~isfield(options,'phiFactor')
                options.phiFactor = defaults.phiFactor;
            end
            if ~isfield(options,'ILambdaDecay')
                options.ILambdaDecay = defaults.ILambdaDecay;
            end
            if ~isfield(options,'penaltyOnly')
                options.penaltyOnly = defaults.penaltyOnly;
            end
            if ~isfield(options,'BoxQPActive')
                options.BoxQPActive = defaults.BoxQPActive;
            end
            if ~isfield(options,'BoxQPBounds')
                options.BoxQPBounds = defaults.BoxQPBounds;
            end
            if ~isfield(options,'SigmaScale')
                options.SigmaScale = defaults.SigmaScale;
            end
            if ~isfield(options,'receedingHorizon')
                options.receedingHorizon = defaults.receedingHorizon;
            end
            if ~isfield(options,'debug')
                options.debug = defaults.debug;
            end
            if ~isfield(options,'derivative_method')
                options.derivative_method = defaults.derivative_method;
            end
            if ~isfield(options,'integration_method')
                options.integration_method = defaults.integration_method;
            end
            
            obj = obj@IndirectTrajectoryOptimization(plant,N,duration,options);
            
            %Figure out which integrator to use
            switch obj.options.integration_method
                case CUDPTrajectoryOptimization.EULER
                    obj.integrator = @obj.euler;
                    obj.integrator_b = @obj.euler_b;
                case CUDPTrajectoryOptimization.MIDPOINT
                    obj.integrator = @obj.midpoint;
                    obj.integrator_b = @obj.midpoint_b;
                case CUDPTrajectoryOptimization.RK3
                    obj.integrator = @obj.rk3;
                    obj.integrator_b = @obj.rk3_b;
                case CUDPTrajectoryOptimization.RK4
                    obj.integrator = @obj.rk4;
                    obj.integrator_b = @obj.rk4_b;
                otherwise
                    error('Drake:CUDPTrajectoryOptimization:InvalidArgument','Unknown integration method');
            end
            
            % initialize to no constraints
            obj.numXinConstraints = 0;
            obj.numXeqConstraints = 0;
            obj.numUinConstraints = 0;
            obj.numUeqConstraints = 0;
            
            % initialize to no callback
            obj.planVisualizer = struct;
            
            % initialize to default receeding multiple function
            obj.computeReceedingMultiple = @obj.computeReceedingMultipleDefault;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Add Function Handles Wrappers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = addPlanVisualizer(obj,callbackFunc,indicies)
            % Add a callback function for a plan visualizer
            obj.planVisualizer.callback = callbackFunc;
            obj.planVisualizer.indicies = indicies;
        end
        function obj = addReceedingMultipleFunction(obj,callbackFunc)
            % Add a callback function for a reeceeding multiple
            obj.computeReceedingMultiple = callbackFunc;
        end
        function obj = addInputConstraint(obj,constraint,Nc,EqFlag)
            % Add constraint (or composite constraint) that is a function of the input
            % @param constraint  a Constraint or CompositeConstraint
            % @param Nc  number of constraints (size of vector returned by
            % constraint function
            if (nargin < 4)
                EqFlag = 0;
            end
            if (EqFlag)
                obj.ueqConstraint = constraint;
                obj.numUeqConstraints = Nc;
            else
                obj.uinConstraint = constraint;
                obj.numUinConstraints = Nc;
            end
        end
        function obj = addStateConstraint(obj,constraint,Nc,EqFlag)
            % Add constraint (or composite constraint) that is a function of the state
            % @param constraint  a CompositeConstraint
            % @param Nc  number of constraints (size of vector returned by
            % constraint function
            if (nargin < 4)
                EqFlag = 0;
            end
            if (EqFlag)
                obj.xeqConstraint = constraint;
                obj.numXeqConstraints = Nc;
            else
                obj.xinConstraint = constraint;
                obj.numXinConstraints = Nc;
            end
        end
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main solve algorithm loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [xtraj,utraj,data] = solveTraj(obj,t_init,traj_init)
            
            n = obj.plant.getNumContStates();
            m = obj.plant.getNumInputs();
            N = obj.getN();
            
            %Allocate matrices
            [mu,l,phi,Id,data] = obj.allocateMatricies();
            
            %Set Initial Trajs and times
            x = zeros(n, N);
            u = zeros(m, N-1);
            x(:,1) = traj_init.x.eval(0);
            if isfield(traj_init,'u')
                u = traj_init.u.eval(traj_init.u.getBreaks());
            end
            t = zeros(N, 1);
            h = (t_init(1)/(N-1))*ones(N-1, 1);
            
            % Set initial feedback to zero
            du = zeros(m,N-1);
            K = zeros(m,n,N-1);
            
            % set initial parameters
            alpha = 1;
            rho = 1;
            drho = 1;
            g_norm = inf;
            deltaJ = inf;
            deltaU = inf;
            maxMu = obj.options.muInit;
            exitPass = 0;
            
            %Initial forward pass to evaluate derivatives of dynamics, cost, and constraints wheere needed
            if obj.options.derivative_method == CUDPTrajectoryOptimization.iLQR
                derivFlag = 1;
            elseif obj.options.derivative_method == CUDPTrajectoryOptimization.UDP
                derivFlag = 0;
            else
                error('Drake:CUDPTrajectoryOptimization:InvalidArgument','Unknown derivative method');
            end
            [x,u,J,Jl,Jmu,A,B,Q,q,Cin,cin,Ceq,ceq,Din,din,Deq,deq,I,maxc,maxd] = obj.forwardPass(t,x,u,h,du,K,...
                                                                                                 alpha,mu,l,derivFlag);
            
            % start major iterations
            majorIter = 0;
            while 1
                %For the exit pass simply continue else
                if ~exitPass
                    %Check iter to see if we should stop for iters
                    if majorIter < obj.options.maxMajorIter
                        majorIter = majorIter+1;
                    else
                        if obj.options.debug
                            disp('Stoping for max majorIter');
                        end
                        break;
                    end
                    % restart minor iter counter
                    iter = 0;
                end
                % do minor iteration
                while 1
                    %Check iter to see if we should stop for iters
                    if iter < obj.options.maxMinorIter
                        iter = iter+1;
                    else
                        if obj.options.debug
                            disp('Stoping for maxMinorIter');
                        end
                        break;
                    end

                    %Do backward pass to calculate du and K
                    done = 0;
                    while ~done
                        if obj.options.derivative_method == CUDPTrajectoryOptimization.iLQR
                            [dunew,Knew,Pnew,pnew,cholFail] = obj.backwardPassiLQR(t,x,u,A,B,Q,q,Cin,cin,Ceq,ceq,...
                                                                                   Din,din,Deq,deq,I,mu,l,rho);
                        elseif obj.options.derivative_method == CUDPTrajectoryOptimization.UDP
                            [dunew,Knew,Pnew,pnew,cholFail] = obj.backwardPassUDP(t,x,u,h,Q,q,Cin,cin,Ceq,ceq,...
                                                                                  Din,din,Deq,deq,I,mu,l,rho);
                        else
                            error('Drake:CUDPTrajectoryOptimization:InvalidArgument','Unknown derivative method');
                        end
                        if cholFail
                            %Increase rho
                            drho = max(drho*obj.options.rhoFactor, obj.options.rhoFactor);
                            rho = max(rho*drho, obj.options.rhoMin);
                            if rho > obj.options.rhoMax
                                break;
                            end
                        else
                            K = Knew;
                            du = dunew;
                            P = Pnew; % broken out for debug purposes at any time
                            p = pnew; % broken out for debug purposes at any time
                            done = 1;
                        end
                    end

                    %Do forward pass with line search to find new J, x, and u
                    alpha = 1;
                    done = 0;
                    while ~done
                        [xnew,unew,Jnew,Jlnew,Jmunew,Anew,Bnew,Qnew,qnew,Cinnew,cinnew,Ceqnew,ceqnew,...
                        Dinnew,dinnew,Deqnew,deqnew,Inew,maxcnew,maxdnew] = obj.forwardPass(t,x,u,h,du,K,...
                                                                                            alpha,mu,l,derivFlag);
                        % compute total costs for comparison
                        Jt = J + Jl + Jmu;
                        Jtnew = Jnew + Jlnew + Jmunew;
                        deltaJ = Jt - Jtnew;
                        deltaU = max(max(abs(u - unew)));
                        
                        % test for acceptance of the new simulation
                        if deltaJ > 0
                            done = 1;
                            g_norm = max(max(abs(du)./(abs(u)+1),[],1));
                            % and save all new vars
                            [x,u,J,Jl,Jmu,maxc,maxd,A,B,Q,q,Cin,cin,Ceq,ceq,Din,din,Deq,deq,I] = ...
                                obj.saveNewVars(xnew,unew,Jnew,Jlnew,Jmunew,Anew,Bnew,Qnew,qnew,...
                                Cinnew,cinnew,Ceqnew,ceqnew,Dinnew,dinnew,Deqnew,deqnew,Inew,maxcnew,maxdnew);
                            %Try to decrease rho if not exit pass
                            drho = min(drho/obj.options.rhoFactor, 1/obj.options.rhoFactor);
                            rho = rho*drho*(rho > obj.options.rhoMin);
                        elseif alpha > obj.options.alphaMin
                            % Failed - decrease alpha first
                            alpha = alpha / obj.options.alphaFactor;
                        else
                            % Failed - increase rho as last resort
                            deltaJ = obj.options.tolFun;
                            drho = max(drho*obj.options.rhoFactor, obj.options.rhoFactor);
                            rho = max(rho*drho, obj.options.rhoMin);
                            done = 1;
                        end
                    end
                    
                    % Save data for this trajectory and print debug
                    maxCon = max(maxc,maxd);
                    if obj.options.debug
                        [maxMu] = max([mu.xin(:);mu.xeq(:);mu.uin(:);mu.ueq(:)]);
                        totalC = sum(abs(ceq(:)))/max(size(ceq)) + sum(cin(:).*(cin(:)>0))/max(size(cin));
                        totalD = sum(abs(deq(:)))/max(size(deq)) + sum(din(:).*(din(:)>0))/max(size(din));
                        [data] = obj.saveDataDebugMinor(data,x,iter,J,Jl,Jmu,maxMu,maxc,maxd,totalC,totalD,deltaJ,deltaU,alpha,N);
                    end
                    
                    %Check deltaJ to see if we should stop
                    if (deltaJ <= (obj.options.tolFunExit*exitPass + obj.options.tolFun*(~exitPass)))
                        if obj.options.debug
                            disp('Stopping for tolFun');
                        end
                        break;
                    end

                    %Check gradient to see if we should stop
                    if (g_norm <= (obj.options.tolGradExit*exitPass + obj.options.tolGrad*(~exitPass)))
                        if obj.options.debug
                            disp('Stopping for tolFun');
                        end
                        break;
                    end

                    %Check regularization parameter to see if we should stop
                    if rho > obj.options.rhoMax
                        if obj.options.debug
                            disp('Stopping for rho');
                        end
                        break;
                    end
                    
                    % If applicable send to plan visualizer
                    if isfield(obj.planVisualizer,'callback')
                        obj.planVisualizer.callback(x(obj.planVisualizer.indicies,:));
                    end
                end
                
                % save data and print debug
                if obj.options.debug
                    [data] = obj.saveDataDebugMajor(data,mu.xin,mu.xeq,mu.uin,mu.ueq,maxMu,J,Jl,Jmu,maxc,maxd,maxCon,iter,majorIter,exitPass);
                end
                
                % If applicable send to plan visualizer
                if isfield(obj.planVisualizer,'callback')
                    obj.planVisualizer.callback(x(obj.planVisualizer.indicies,:));
                end
                
                % see if we should exit
                if obj.options.receedingHorizon
                    readyToExit = obj.receedingHorizonCheck(cin,ceq,din,deq,xinTol,xeqTol,uinTol,ueqTol);
                else
                    readyToExit = maxCon < obj.options.tolConstr;
                end
                % make sure we do an exit pass to hit exit tolerance
                if readyToExit
                    if exitPass
                        break;
                    else
                        exitPass = 1;
                    end
                else
                    % else prepare for next major iteration
                    exitPass = 0;
                    [mu,l,phi,Id,Jl,Jmu] = muLambdaUpdate(obj,mu,l,phi,I,Id,cin,ceq,din,deq);
                end
            end
            % clean up and return
            xtraj = obj.reconstructStateTrajectory(x);
            if nargout>1
                utraj = obj.reconstructInputTrajectory(u);
            end
        end
    end
    
    methods (Access = private)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAIN FORWARD BACKWARD PASSES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Forward pass with simulation
        function [xnew,unew,Jnew,Jlnew,Jmunew,A,B,Q,q,Cin,cin,Ceq,ceq,Din,din,Deq,deq,I,maxc,maxd] = forwardPass(obj,t,x,u,h,du,K,alpha,mu,l,derivFlag)
            n = size(x,1);
            m = size(u,1);
            N = size(x,2);
            % allocate matricies
            unew = zeros(m,N-1);
            xnew = zeros(n,N);
            Q = zeros(n+m,n+m,N);
            q = zeros(n+m,N);
            A = zeros(n,n,N-1);
            B = zeros(n,m,N-1);
            % set initial point
            xnew(:,1) = x(:,1);
            Jnew = 0;
            
            % then loop thorugh and integrate forward
            for k = 1:N-1
                unew(:,k) = u(:,k) - alpha*du(:,k) - K(:,:,k)*(xnew(:,k)-x(:,k));
                if derivFlag
                    [xnew(:,k+1),A(:,:,k),B(:,:,k)] = obj.integrator(xnew(:,k),unew(:,k),h(k));
                else
                    xnew(:,k+1) = obj.integrator(xnew(:,k),unew(:,k),h(k));
                end
                if (any(isnan(xnew(:,k+1))) || any(isinf(xnew(:,k+1))))
                    Jnew = inf;
                    [Cin,cin,Ceq,ceq,Din,din,Deq,deq] = obj.allocateCnstrMatricies(m,n,N);
                    I = struct;
                    Jlnew = inf;
                    Jmunew = inf;
                    maxc = inf;
                    maxd = inf;
                    return;
                end
                [cost, q(:,k), Q(:,:,k)] = obj.running_cost(h(k),xnew(:,k),unew(:,k));
                Jnew = Jnew + cost;
            end
            [cost, q(1:n,N), Q(1:n,1:n,N)] = obj.final_cost(t(N), xnew(:,N));
            Jnew = Jnew + cost;
            
            % and compute constraints and their costs
            [Cin,cin,Ceq,ceq,Din,din,Deq,deq,I,Jlnew,Jmunew,maxc,maxd] = obj.computeConstraints(xnew,unew,mu,l);
        end
        
        % UDP LQR update back in time
        function [du,K,P,p,cholFail] = backwardPassUDP(obj,t,x,u,h,Q,q,Cin,cin,Ceq,ceq,Din,din,Deq,deq,I,mu,l,rho)
            
            %Set up backwards LQR pass
            n = size(x,1);
            m = size(u,1);
            N = size(x,2);
            P = zeros(n,n,N);
            p = zeros(n, N);
            K = zeros(m,n,N-1);
            du = zeros(m,N-1);
            
            % Last step cost
            P(:,:,N) = Q(1:n,1:n,N);
            p(:,N) = q(1:n,N);
            % Add constraints
            P(:,:,N) = P(:,:,N) + Cin(:,:,N)'*diag(mu.xin(:,N))*diag(I.xin(:,N))*Cin(:,:,N) + ...
                                  Ceq(:,:,N)'*diag(mu.xeq(:,N))*Ceq(:,:,N);
            p(:,N) = p(:,N) + Cin(:,:,N)'*diag(mu.xin(:,N))*diag(I.xin(:,N))*cin(:,N) + ...
                              Ceq(:,:,N)'*diag(mu.xeq(:,N))*ceq(:,N) + ...
                              Cin(:,:,N)'*l.xin(:,N) + Ceq(:,:,N)'*l.xeq(:,N);
            
            for k = (N-1):-1:1
                % Compute running cost and penalty + lambda terms for this time step
                currH = Q(:,:,k) + blkdiag(Cin(:,:,k)'*diag(mu.xin(:,k))*diag(I.xin(:,k))*Cin(:,:,k) + ...
                                           Ceq(:,:,k)'*diag(mu.xeq(:,k))*Ceq(:,:,k), ...
                                           Din(:,:,k)'*diag(mu.uin(:,k))*diag(I.uin(:,k))*Din(:,:,k) + ...
                                           Deq(:,:,k)'*diag(mu.ueq(:,k))*Deq(:,:,k));
                currg = q(:,k) + [Cin(:,:,k)'*diag(mu.xin(:,k))*diag(I.xin(:,k))*cin(:,k) + ...
                                  Ceq(:,:,k)'*diag(mu.xeq(:,k))*ceq(:,k); ...
                                  Din(:,:,k)'*diag(mu.uin(:,k))*diag(I.uin(:,k))*din(:,k) + ...
                                  Deq(:,:,k)'*diag(mu.ueq(:,k))*deq(:,k)] + ...
                                  [Cin(:,:,k)'*l.xin(:,k) + Ceq(:,:,k)'*l.xeq(:,k); ...
                                   Din(:,:,k)'*l.uin(:,k) + Deq(:,:,k)'*l.ueq(:,k)];
                
                %Generate sigma points from Pxx(k+1) and Huu(k)
                [S,cholFail] = chol(blkdiag(inv(P(1:n,1:n,k+1)), inv(currH((n+1):end,(n+1):end))), 'lower');
                if cholFail
                    return;
                end
                S = obj.options.SigmaScale*S;
                Sig = [S -S];
                for j = 1:(2*(n+m))
                    Sig(:,j) = Sig(:,j) + [x(:,k+1); u(:,k)];
                end

                %Project px(k+1) onto sigma points
                g = zeros(2*(n+m),1);
                for j = 1:(n+m)
                    g(j) = p(1:n,k+1)'*S(1:n,j);
                    g(j+n+m) = -g(j);
                end

                %Propagate sigma points through backwards dynamics
                for j = 1:(2*(n+m))
                    Sig(1:n,j) = obj.integrator_b(Sig(1:n,j),Sig(n+1:end,j),h(k));
                end

                %Calculate [gx; gu] from sigma points
                D = zeros(n+m,n+m);
                df = zeros(n+m,1);
                for j = 1:(n+m)
                    D(j,:) = (Sig(:,j)-Sig(:,n+m+j))';
                    df(j) = g(j) - g(n+m+j);
                end
                gt = D\df;
                gx = gt(1:n) + currg(1:n); %add on one-step cost
                gu = gt((n+1):end) + currg(n+1:end); %add on one-step cost

                %Calculate Hessian w.r.t. [x_k; u_k] from sigma points
                eta = zeros(n+m,1);
                for j = 1:(2*(n+m))
                    eta = eta + (1/(2*(n+m)))*Sig(:,j);
                end
                M = zeros(n+m);
                for j = 1:(2*(n+m))
                    M = M + (0.5/(obj.options.SigmaScale^2))*(Sig(:,j)-eta)*(Sig(:,j)-eta)';
                end
                H = M^-1;
                H(1:n,1:n) = H(1:n,1:n) + currH(1:n,1:n); %add in one-step state cost for this timestep
                H(n+1:end,1:n) = H(n+1:end,1:n) + currH(n+1:end,1:n);
                H(1:n,n+1:end) = H(1:n,n+1:end) + currH(1:n,n+1:end);

                Hxx = H(1:n,1:n);
                Huu = H(n+1:end,n+1:end);
                Hux = H(n+1:end,1:n);
                Hxu = H(1:n,n+1:end);
                HuuReg = Huu + rho*eye(m); % add in regularization term

                if obj.options.BoxQPActive
                    % Solve QP to find du and K
                    lb = obj.options.BoxQPBounds(1) - u(:,k);
                    ub = obj.options.BoxQPBounds(2) - u(:,k);
                    [du(:,k),Pfree,Lfree,cholFail] = boxQP(HuuReg,gu,lb,ub,du(:,k));
                    if cholFail
                        break;
                    end    
                    K(:,:,k) = zeros(m,n);
                    Kfree = Lfree'\(Lfree\(Pfree*Hux));
                    Kfree(~isfinite(Kfree)) = 0;
                    K(:,:,k) = Pfree'*Kfree;
                    du(:,k) = -1*du(:,k);
                else
                    % Check for positive definiteness
                    [Suu,cholFail] = chol(HuuReg);
                    if cholFail
                        break;
                    end
                    % Calculate du and K
                    du(:,k) = Suu\(Suu'\gu);
                    K(:,:,k) = Suu\(Suu'\Hux);
                end
                
                %Calculate new cost-to-go function
                p(:,k) = -Hxu*du(:,k) + K(:,:,k)'*Huu*du(:,k) + gx - K(:,:,k)'*gu;
                P(:,:,k) = Hxx + K(:,:,k)'*Huu*K(:,:,k) - Hxu*K(:,:,k) - K(:,:,k)'*Hux;
            end
        end

        % LQR update back in time
        function [du,K,P,p,cholFail] = backwardPassiLQR(obj,t,x,u,A,B,Q,q,Cin,cin,Ceq,ceq,Din,din,Deq,deq,I,mu,l,rho)
            
            %Set up backwards LQR pass
            n = size(x,1);
            m = size(u,1);
            N = size(x,2);
            P = zeros(n,n,N);
            p = zeros(n, N);
            K = zeros(m,n,N-1);
            du = zeros(m,N-1);
            
            % Last step cost
            P(:,:,N) = Q(1:n,1:n,N);
            p(:,N) = q(1:n,N);
            % Add constraints
            P(:,:,N) = P(:,:,N) + Cin(:,:,N)'*diag(mu.xin(:,N))*diag(I.xin(:,N))*Cin(:,:,N) + ...
                                  Ceq(:,:,N)'*diag(mu.xeq(:,N))*Ceq(:,:,N);
            p(:,N) = p(:,N) + Cin(:,:,N)'*diag(mu.xin(:,N))*diag(I.xin(:,N))*cin(:,N) + ...
                              Ceq(:,:,N)'*diag(mu.xeq(:,N))*ceq(:,N) + ...
                              Cin(:,:,N)'*l.xin(:,N) + Ceq(:,:,N)'*l.xeq(:,N);
            
            for k = (N-1):-1:1
                % Propagate cost-to-go function through linearized dynamics
                H = [A(:,:,k) B(:,:,k)]'*P(:,:,k+1)*[A(:,:,k) B(:,:,k)];
                g = [A(:,:,k) B(:,:,k)]'*p(:,k+1);
                
                % Add on running cost for this time step
                H = H + Q(:,:,k);
                g = g + q(:,k);
                
                % Add on penalty + lambda terms for this time step
                H = H + blkdiag(Cin(:,:,k)'*diag(mu.xin(:,k))*diag(I.xin(:,k))*Cin(:,:,k) + ... 
                                Ceq(:,:,k)'*diag(mu.xeq(:,k))*Ceq(:,:,k), ...
                                Din(:,:,k)'*diag(mu.uin(:,k))*diag(I.uin(:,k))*Din(:,:,k) + ...
                                Deq(:,:,k)'*diag(mu.ueq(:,k))*Deq(:,:,k));
                g = g + [Cin(:,:,k)'*diag(mu.xin(:,k))*diag(I.xin(:,k))*cin(:,k) + ...
                         Ceq(:,:,k)'*diag(mu.xeq(:,k))*ceq(:,k); ...
                         Din(:,:,k)'*diag(mu.uin(:,k))*diag(I.uin(:,k))*din(:,k) + ...
                         Deq(:,:,k)'*diag(mu.ueq(:,k))*deq(:,k)] + ...
                         [Cin(:,:,k)'*l.xin(:,k) + Ceq(:,:,k)'*l.xeq(:,k); ...
                          Din(:,:,k)'*l.uin(:,k) + Deq(:,:,k)'*l.ueq(:,k)];
                
                % Break cost to go matrix/vector into convenient blocks
                Hxx = H(1:n,1:n);
                Huu = H((n+1):end,(n+1):end);
                Hxu = H(1:n,(n+1):end);
                Hux = H((n+1):end,1:n);
                gx = g(1:n);
                gu = g((n+1):end);
                
                % Add regularization term
                HuuReg = Huu + rho*eye(m);
                
                if obj.options.BoxQPActive
                    % Solve QP to find du and K
                    lb = obj.options.BoxQPBounds(1) - u(:,k);
                    ub = obj.options.BoxQPBounds(2) - u(:,k);
                    [du(:,k),Pfree,Lfree,cholFail] = boxQP(HuuReg,gu,lb,ub,du(:,k));
                    if cholFail
                        break;
                    end    
                    K(:,:,k) = zeros(m,n);
                    Kfree = Lfree'\(Lfree\(Pfree*Hux));
                    Kfree(~isfinite(Kfree)) = 0;
                    K(:,:,k) = Pfree'*Kfree;
                    du(:,k) = -1*du(:,k);
                else
                    % Check for positive definiteness
                    [Suu,cholFail] = chol(HuuReg);
                    if cholFail
                        break;
                    end
                    % Calculate du and K
                    du(:,k) = Suu\(Suu'\gu);
                    K(:,:,k) = Suu\(Suu'\Hux);
                end
                
                %Calculate new cost-to-go function
                p(:,k) = -Hxu*du(:,k) + K(:,:,k)'*Huu*du(:,k) + gx - K(:,:,k)'*gu;
                P(:,:,k) = Hxx + K(:,:,k)'*Huu*K(:,:,k) - Hxu*K(:,:,k) - K(:,:,k)'*Hux;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MULTIPLIER UPDATE FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % update mu and lambda for a single variable type (xin vs ueq) vectorized
        function [mu,l,phi,Jl,Jmu,Id] = updateMuLambdaSingle(obj,mu,l,phi,Jl,Jmu,cd,receedingFactor,I,Id)
            % if I,Id passed in we are in an inequality situation
            if nargin > 8
                eqFlag = 0;
            else
                eqFlag = 1;
                Id = 0;
            end
            % resize the receeding factor for vector ops
            receedingFactor = repmat(receedingFactor,size(cd,1),1);

            % if penalty only we only update mu and Jmu
            if obj.options.penaltyOnly
                muFlags = double(mu < obj.options.muMax);
            else
                % else we need to figure out if we are below phi for lambdas or just do
                % the update on the mus by first computing flags to indicate which to update
                % note that on equality it is absolute not signed distance that matters
                error = (cd./receedingFactor);
                if eqFlag
                    error = abs(error);
                end
                lambdaFlags = (error <= phi);
                muFlags = (mu < obj.options.muMax).*~lambdaFlags;
                % update lambda first because we need the old mu to do so
                if eqFlag
                    l = l + mu.*lambdaFlags.*cd;
                % if inequality take into account the I and Id
                else
                    Id = max(I, Id - 1/obj.options.ILambdaDecay);
                    l = l + mu.*lambdaFlags.*Id.*cd;
                    l = max(0, l);
                end
                % and update Jl
                Jl = Jl + sum(sum(l.*cd,1),2);
                % and update the phierance matrix
                phi = max(1/obj.options.phiFactor.*lambdaFlags.*phi + ...
                          ~lambdaFlags.*phi, obj.options.tolConstr);
            end
            % in either case then update mu and Jmu
            mu = obj.options.muFactor.*muFlags.*mu + ~muFlags.*mu;
            if eqFlag
                Jmu = Jmu + sum(sum(cd.*mu.*cd,1),2);
            % if inequality take into account the I again
            else
                Jmu = Jmu + sum(sum(cd.*mu.*I.*cd,1),2);
            end
        end
        
        % update mu and lambda based on constraint violations
        function [mu,l,phi,Id,Jl,Jmu] = muLambdaUpdate(obj,mu,l,phi,I,Id,cin,ceq,din,deq)
            N = obj.getN();
            Jl = 0;
            Jmu = 0;
            
            % test each constraint on each timestep in vectorized form for speed
            % first compute the reeceding multiple if that is turned on
            if obj.options.receedingHorizon
                receedingFactor = arrayfun(@obj.computeReceedingMultiple,1:N);
            else
                receedingFactor = ones(1,N);
            end
            % then update for each type of constraint (if they exist)
            if obj.numXinConstraints
                [mu.xin,l.xin,phi.xin,Jl,Jmu,Id.xin] = obj.updateMuLambdaSingle(...
                    mu.xin,l.xin,phi.xin,Jl,Jmu,cin,receedingFactor,I.xin,Id.xin);
            end
            if obj.numXeqConstraints
                [mu.xeq,l.xeq,phi.xeq,Jl,Jmu,~] = obj.updateMuLambdaSingle(...
                    mu.xeq,l.xeq,phi.xeq,Jl,Jmu,ceq,receedingFactor);
            end
            if obj.numUinConstraints
                [mu.uin,l.uin,phi.uin,Jl,Jmu,Id.uin] = obj.updateMuLambdaSingle(...
                    mu.uin,l.uin,phi.uin,Jl,Jmu,din,receedingFactor(:,1:end-1),I.uin,Id.uin);
            end
            if obj.numUeqConstraints
                [mu.ueq,l.ueq,phi.ueq,Jl,Jmu,~] = obj.updateMuLambdaSingle(...
                    mu.ueq,l.ueq,phi.ueq,Jl,Jmu,deq,receedingFactor(:,1:end-1));
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Debug and Other (non Integrators) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Compute constraints
        function [Cin,cin,Ceq,ceq,Din,din,Deq,deq,I,Jlnew,Jmunew,maxc,maxd] = computeConstraints(obj,xnew,unew,mu,l)
            n = size(xnew,1);
            m = size(unew,1);
            N = size(xnew,2);
            % zero out constraint indicators and costs
            I.xin = zeros(obj.numXinConstraints, N);
            I.uin = zeros(obj.numUinConstraints, N-1);
            Jlnew = 0;
            Jmunew = 0;
            maxc = 0;
            maxd = 0;
            % allocate matricies
            [Cin,cin,Ceq,ceq,Din,din,Deq,deq] = obj.allocateCnstrMatricies(m,n,N);
            for k = 1:N
                % if constarints include them
                if (obj.numXinConstraints)
                    [cin(:,k), Cin(:,:,k)] = obj.xinConstraint(k,xnew(:,k));
                    I.xin(:,k) = (cin(:,k) >= 0) | (l.xin(:,k) > 0);
                    if max(cin(:,k)) > maxc
                        maxc = max(cin(:,k));
                    end
                    Jlnew = Jlnew + l.xin(:,k)'*cin(:,k);
                    Jmunew = Jmunew + 0.5*cin(:,k)'*diag(mu.xin(:,k))*diag(I.xin(:,k))*cin(:,k);
                end
                if (obj.numXeqConstraints)
                    [ceq(:,k), Ceq(:,:,k)] = obj.xeqConstraint(k,xnew(:,k));
                    if max(abs(ceq(:,k))) > maxc
                        maxc = max(abs(ceq(:,k)));
                    end
                    Jlnew = Jlnew + l.xeq(:,k)'*ceq(:,k);
                    Jmunew = Jmunew + 0.5*ceq(:,k)'*diag(mu.xeq(:,k))*ceq(:,k);
                end
                if (obj.numUinConstraints && k < N)
                    [din(:,k), Din(:,:,k)] = obj.uinConstraint(k,unew(:,k));
                    I.uin(:,k) = (din(:,k) >= 0) | (l.uin(:,k) > 0);
                    if max(din(:,k)) > maxd
                        maxd = max(din(:,k));
                    end
                    Jlnew = Jlnew + l.uin(:,k)'*din(:,k);
                    Jmunew = Jmunew + 0.5*din(:,k)'*diag(mu.uin(:,k))*diag(I.uin(:,k))*din(:,k);
                end
                if (obj.numUeqConstraints && k < N)
                    [deq(:,k), Deq(:,:,k)] = obj.ueqConstraint(k,unew(:,k));
                    if max(abs(deq(:,k))) > maxd
                        maxd = max(abs(deq(:,k)));
                    end
                    Jlnew = Jlnew + l.ueq(:,k)'*deq(:,k);
                    Jmunew = Jmunew + 0.5*deq(:,k)'*diag(mu.ueq(:,k))*deq(:,k);
                end
            end
        end
        
        % Debug Functions
        function [data] = saveDataDebugMinor(obj,data,x,iter,J,Jl,Jmu,maxMu,maxc,maxd,totalC,totalD,deltaJ,deltaU,alpha,N)
            data.Jinner(end+1) = J;
            data.Jtinner(end+1) = J + Jl + Jmu;
            data.cinner(end+1) = maxc;
            data.dinner(end+1) = maxd;
            data.ainner(end+1) = alpha;
            data.dJinner(end+1) = deltaJ;
            data.totalC(end+1) = totalC;
            data.totalD(end+1) = totalD;
            disp(['iter [',num2str(iter),'] xN[',num2str(x(1,N)),' ',num2str(x(2,N)),'] log(max(mu))[',num2str(log10(maxMu)),'] deltaU[',num2str(deltaU),'] dJ[',num2str(deltaJ),'] maxc/d[',num2str(maxc),' ',num2str(maxd),']']);
        end
        function [data] = saveDataDebugMajor(obj,data,muXin,muXeq,muUin,muUeq,maxMu,J,Jl,Jmu,maxc,maxd,maxCon,iter,majorIter,exitPass)
            if (obj.numXinConstraints)
                data.muXin(end+1) = max(muXin(:));
            end
            if (obj.numXeqConstraints)
                data.muXeq(end+1) = max(muXeq(:));
            end
            if (obj.numUinConstraints)
                data.muUin(end+1) = max(muUin(:));
            end
            if (obj.numUeqConstraints)
                data.muUeq(end+1) = max(muUeq(:));
            end
            data.majorIter(end+1) = majorIter;
            data.iter(end+1) = iter;
            data.J(end+1) = J;
            data.Jt(end+1) = J + Jl + Jmu;
            data.maxCon(end+1) = maxCon;
            data.maxc(end+1) = maxc;
            data.maxd(end+1) = maxd;
            data.exitPass(end+1) = exitPass;
            disp('-------------------------------');
            disp(['iter [',num2str(majorIter),'] log(max(mu))[',num2str(log10(maxMu)),'] J[',num2str(J + Jl + Jmu),'] maxCon[',num2str(max(maxc,maxd)),']']);
            disp('-------------------------------');
        end
        
        % Receeding Horizon Functions
        function [factor] = computeReceedingMultipleDefault(obj,k)
            alpha = 0.5;
            if k > 1
                factor = exp(alpha*k);
            else
                factor = 1;
            end
        end
        function [readyToExit] = receedingHorizonCheck(obj,cin,ceq,din,deq,lxin,lxeq,luin,lueq)
            readyToExit = 1;
            N = obj.getN();
            for k = 1:N
                factor = obj.computeReceedingMultiple(k);
                if (any(cin(:,k)/factor > lxin(:,k)) || any(ceq(:,k)/factor > lxeq(:,k)))
                    readyToExit = 0;
                    break;
                end
                if (k < N)
                    if (any(din(:,k)/factor > luin(:,k)) || any(deq(:,k)/factor > lueq(:,k)))
                        readyToExit = 0;
                        break;
                    end
                end
            end
        end
        
        % Initialize empty variables for mulitpliers
        function [mu,l,phi,Id,data] = allocateMatricies(obj)
            N = obj.getN();
            
            Id.xin = zeros(obj.numXinConstraints, N);
            Id.uin = zeros(obj.numUinConstraints, N-1);
            l.xin = zeros(obj.numXinConstraints, N);
            l.xeq = zeros(obj.numXeqConstraints, N);
            l.uin = zeros(obj.numUinConstraints, N-1);
            l.ueq = zeros(obj.numUeqConstraints, N-1);
            mu.xin = obj.options.muInit * ones(obj.numXinConstraints, N);
            mu.xeq = obj.options.muInit * ones(obj.numXeqConstraints, N);
            mu.uin = obj.options.muInit * ones(obj.numUinConstraints, N-1);
            mu.ueq = obj.options.muInit * ones(obj.numUeqConstraints, N-1);
            phi.xin = obj.options.tolLambdaUpdate * ones(obj.numXinConstraints, N);
            phi.xeq = obj.options.tolLambdaUpdate * ones(obj.numXeqConstraints, N);
            phi.uin = obj.options.tolLambdaUpdate * ones(obj.numUinConstraints, N-1);
            phi.ueq = obj.options.tolLambdaUpdate * ones(obj.numUeqConstraints, N-1);
            data.Jinner = [];
            data.Jtinner = [];
            data.cinner = [];
            data.dinner = [];
            data.ainner = [];
            data.dJinner = [];
            data.totalC = [];
            data.totalD = [];
            data.muXin = [];
            data.muXeq = [];
            data.muUin = [];
            data.muUeq = [];
            data.majorIter = [];
            data.iter = [];
            data.J = [];
            data.Jt = [];
            data.maxCon = [];
            data.maxc = [];
            data.maxd = [];
            data.exitPass = [];
        end
        % And for constraints
        function [Cin,cin,Ceq,ceq,Din,din,Deq,deq] = allocateCnstrMatricies(obj,m,n,N)
            Cin = zeros(obj.numXinConstraints, n, N);
            cin = zeros(obj.numXinConstraints, N);
            Ceq = zeros(obj.numXeqConstraints, n, N);
            ceq = zeros(obj.numXeqConstraints, N);
            Din = zeros(obj.numUinConstraints, m, N-1);
            din = zeros(obj.numUinConstraints, N-1);
            Deq = zeros(obj.numUeqConstraints, m, N-1);
            deq = zeros(obj.numUeqConstraints, N-1);
        end
        
        % save new variables over to old names by just assignign outputs to inputs
        function [x,u,J,Jl,Jmu,maxc,maxd,A,B,Q,q,Cin,cin,Ceq,ceq,Din,din,Deq,deq,I] = ...
                 saveNewVars(obj,x,u,J,Jl,Jmu,A,B,Q,q,Cin,cin,Ceq,ceq,Din,din,Deq,deq,I,maxc,maxd)
             return;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INTEGRATORS FOLLOW %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
        % All of these integrators integrate Jacobians as well as states
        function [xkp1,A,B] = euler(obj,xk,uk,hk)
            if nargout == 1
                xdot = obj.plant.dynamics(0,xk,uk);
                xkp1 = xk + hk*xdot;
            else
                [xdot, dxdot] = obj.plant.dynamics(0,xk,uk);
                xkp1 = xk + hk*xdot;
                n = length(xk);
                A = eye(n) + full(dxdot(:,2:(n+1)));
                B = full(dxdot(:,(n+2):end));
            end
        end 
        function [xk,A,B] = euler_b(obj,xkp1,uk,hk)
            if nargout == 1
                xdot = obj.plant.dynamics(0,xkp1,uk);
                xk = xkp1 - hk*xdot;
            else
                [xdot, dxdot] = obj.plant.dynamics(0,xkp1,uk);
                xk = xkp1 - hk*xdot;
                n = length(xkp1);
                A = eye(n) - full(dxdot(:,2:(n+1)));
                B = -1*full(dxdot(:,(n+2):end));
            end
        end
        
        function [xkp1,A,B] = midpoint(obj,xk,uk,hk)
            if nargout == 1
                xdot1 = obj.plant.dynamics(0,xk,uk);
                xdot2 = obj.plant.dynamics(0,xk+.5*hk*xdot1,uk);
                xkp1 = xk + hk*xdot2;
            else
                n = length(xk);
                
                [xdot1, dxdot1] = obj.plant.dynamics(0,xk,uk);
                [xdot2, dxdot2] = obj.plant.dynamics(0,xk+.5*hk*xdot1,uk);
                xkp1 = xk + hk*xdot2;
                
                A1 = eye(n) + (hk/2)*full(dxdot1(:,2:(n+1)));
                B1 = (hk/2)*full(dxdot1(:,(n+2):end));
                
                A2 = eye(n) + (hk/2)*full(dxdot2(:,2:(n+1)));
                B2 = (hk/2)*full(dxdot2(:,(n+2):end));
                
                A = A2*A1;
                B = A2*B1 + B2;
            end
        end 
        function [xk,A,B] = midpoint_b(obj,xkp1,uk,hk)
            if nargout == 1
                xdot1 = obj.plant.dynamics(0,xkp1,uk);
                xdot2 = obj.plant.dynamics(0,xkp1-.5*hk*xdot1,uk);
                xk = xkp1 - hk*xdot2;
            else
                n = length(xkp1);
                
                [xdot1, dxdot1] = obj.plant.dynamics(0,xkp1,uk);
                [xdot2, dxdot2] = obj.plant.dynamics(0,xkp1-.5*hk*xdot1,uk);
                xk = xkp1 - hk*xdot2;
                
                A1 = eye(n) - (hk/2)*full(dxdot1(:,2:(n+1)));
                B1 = -1*(hk/2)*full(dxdot1(:,(n+2):end));
                
                A2 = eye(n) - (hk/2)*full(dxdot2(:,2:(n+1)));
                B2 = -1*(hk/2)*full(dxdot2(:,(n+2):end));
                
                A = A2*A1;
                B = A2*B1 + B2;
            end
        end
        
        function [xkp1,A,B] = rk3(obj,xk,uk,hk)
            if nargout == 1
                xdot1 = obj.plant.dynamics(0,xk,uk);
                xdot2 = obj.plant.dynamics(0,xk+.5*hk*xdot1,uk);
                xdot3 = obj.plant.dynamics(0,xk+.75*hk*xdot2,uk);
                xkp1 = xk + (hk/9)*(2*xdot1 + 3*xdot2 + 4*xdot3);
            else
                n = length(xk);
                
                [xdot1, dxdot1] = obj.plant.dynamics(0,xk,uk);
                [xdot2, dxdot2] = obj.plant.dynamics(0,xk+.5*hk*xdot1,uk);
                [xdot3, dxdot3] = obj.plant.dynamics(0,xk+.75*hk*xdot2,uk);
                xkp1 = xk + (hk/9)*(2*xdot1 + 3*xdot2 + 4*xdot3);
                
                A1 = eye(n) + (2*hk/9)*full(dxdot1(:,2:(n+1)));
                B1 = (2*hk/9)*full(dxdot1(:,(n+2):end));
                
                A2 = eye(n) + (hk/3)*full(dxdot2(:,2:(n+1)));
                B2 = (hk/3)*full(dxdot2(:,(n+2):end));
                
                A3 = eye(n) + (4*hk/9)*full(dxdot3(:,2:(n+1)));
                B3 = (4*hk/9)*full(dxdot3(:,(n+2):end));
                
                A = A3*A2*A1;
                B = A3*A2*B1 + A3*B2 + B3;
            end
        end
        function [xk,A,B] = rk3_b(obj,xkp1,uk,hk)
            if nargout == 1
                xdot1 = obj.plant.dynamics(0,xkp1,uk);
                xdot2 = obj.plant.dynamics(0,xkp1-.5*hk*xdot1,uk);
                xdot3 = obj.plant.dynamics(0,xkp1-.75*hk*xdot2,uk);
                xk = xkp1 - (hk/9)*(2*xdot1 + 3*xdot2 + 4*xdot3);
            else
                n = length(xkp1);
                
                [xdot1, dxdot1] = obj.plant.dynamics(0,xkp1,uk);
                [xdot2, dxdot2] = obj.plant.dynamics(0,xkp1-.5*hk*xdot1,uk);
                [xdot3, dxdot3] = obj.plant.dynamics(0,xkp1-.75*hk*xdot2,uk);
                xk = xkp1 - (hk/9)*(2*xdot1 + 3*xdot2 + 4*xdot3);
                
                A1 = eye(n) - (2*hk/9)*full(dxdot1(:,2:(n+1)));
                B1 = -1*(2*hk/9)*full(dxdot1(:,(n+2):end));
                
                A2 = eye(n) - (hk/3)*full(dxdot2(:,2:(n+1)));
                B2 = -1*(hk/3)*full(dxdot2(:,(n+2):end));
                
                A3 = eye(n) - (4*hk/9)*full(dxdot3(:,2:(n+1)));
                B3 = -1*(4*hk/9)*full(dxdot3(:,(n+2):end));
                
                A = A3*A2*A1;
                B = A3*A2*B1 + A3*B2 + B3;
            end
        end
        
        function [xkp1,A,B] = rk4(obj,xk,uk,hk)
            if nargout == 1
                xdot1 = obj.plant.dynamics(0,xk,uk);
                xdot2 = obj.plant.dynamics(0,xk+.5*hk*xdot1,uk);
                xdot3 = obj.plant.dynamics(0,xk+.5*hk*xdot2,uk);
                xdot4 = obj.plant.dynamics(0,xk+hk*xdot3,uk);
                
                xkp1 = xk + (hk/6)*(xdot1 + 2*xdot2 + 2*xdot3 + xdot4);
            else
                n = length(xk);
                
                [xdot1, dxdot1] = obj.plant.dynamics(0,xk,uk);
                [xdot2, dxdot2] = obj.plant.dynamics(0,xk+.5*hk*xdot1,uk);
                [xdot3, dxdot3] = obj.plant.dynamics(0,xk+.5*hk*xdot2,uk);
                [xdot4, dxdot4] = obj.plant.dynamics(0,xk+hk*xdot3,uk);
                
                xkp1 = xk + (hk/6)*(xdot1 + 2*xdot2 + 2*xdot3 + xdot4);
                
                A1 = eye(n) + (hk/6)*full(dxdot1(:,2:(n+1)));
                B1 = (hk/6)*full(dxdot1(:,(n+2):end));
                
                A2 = eye(n) + (hk/3)*full(dxdot2(:,2:(n+1)));
                B2 = (hk/3)*full(dxdot2(:,(n+2):end));
                
                A3 = eye(n) + (hk/3)*full(dxdot3(:,2:(n+1)));
                B3 = (hk/3)*full(dxdot3(:,(n+2):end));
                
                A4 = eye(n) + (hk/6)*full(dxdot4(:,2:(n+1)));
                B4 = (hk/6)*full(dxdot4(:,(n+2):end));
                
                A = A4*A3*A2*A1;
                B = A4*A3*A2*B1 + A4*A3*B2 + A4*B3 + B4;
            end
        end
        function [xk,A,B] = rk4_b(obj,xkp1,uk,hk)
            if nargout == 1
                xdot1 = obj.plant.dynamics(0,xkp1,uk);
                xdot2 = obj.plant.dynamics(0,xkp1-.5*hk*xdot1,uk);
                xdot3 = obj.plant.dynamics(0,xkp1-.5*hk*xdot2,uk);
                xdot4 = obj.plant.dynamics(0,xkp1-hk*xdot3,uk);
                
                xk = xkp1 - (hk/6)*(xdot1 + 2*xdot2 + 2*xdot3 + xdot4);
            else
                n = length(xkp1);
                
                [xdot1, dxdot1] = obj.plant.dynamics(0,xkp1,uk);
                [xdot2, dxdot2] = obj.plant.dynamics(0,xkp1-.5*hk*xdot1,uk);
                [xdot3, dxdot3] = obj.plant.dynamics(0,xkp1-.5*hk*xdot2,uk);
                [xdot4, dxdot4] = obj.plant.dynamics(0,xkp1-hk*xdot3,uk);
                
                xk = xkp1 - (hk/6)*(xdot1 + 2*xdot2 + 2*xdot3 + xdot4);
                
                A1 = eye(n) - (hk/6)*full(dxdot1(:,2:(n+1)));
                B1 = -1*(hk/6)*full(dxdot1(:,(n+2):end));
                
                A2 = eye(n) - (hk/3)*full(dxdot2(:,2:(n+1)));
                B2 = -1*(hk/3)*full(dxdot2(:,(n+2):end));
                
                A3 = eye(n) - (hk/3)*full(dxdot3(:,2:(n+1)));
                B3 = -1*(hk/3)*full(dxdot3(:,(n+2):end));
                
                A4 = eye(n) - (hk/6)*full(dxdot4(:,2:(n+1)));
                B4 = -1*(hk/6)*full(dxdot4(:,(n+2):end));
                
                A = A4*A3*A2*A1;
                B = A4*A3*A2*B1 + A4*A3*B2 + A4*B3 + B4;
            end
        end
    end 
end
