classdef VariationalTrajectoryOptimization < DirectTrajectoryOptimization
    
    properties
        nC %Number of contact points
        nD %Number of basis vectors in polyhedral friction cone
        nL %Number of contact constraints
        integration_method %Midpoint or Simpson's rule for integration
        angle_inds %Indices of q that are joint angles and should be diffed accounting for wrap-around
        v0_inds
        psi_inds
        eta_inds
        c_inds
        b_inds
        s_inds
    end
    
    properties (Constant)
        MIDPOINT = 2;
        SIMPSON = 3;
    end
    
    methods
        function obj = VariationalTrajectoryOptimization(plant,N,duration,options)
            if nargin < 4
                options=struct();
            end
            
            if ~isfield(options,'active_collision_options')
                options.active_collision_options.terrain_only = true;
            end
            if ~isfield(options,'integration_method')
                options.integration_method = VariationalTrajectoryOptimization.MIDPOINT;
            end
            if ~isfield(options,'multiple_contacts')
                options.multiple_contacts = false;
            end
            if isa(plant,'PlanarRigidBodyManipulator')
                options.twoD = true;
            else
                options.twoD = false;
            end
            
            if ~isfield(options,'s_weight')
                options.s_weight = 10;
            end
            if ~isfield(options,'add_ccost')
                options.add_ccost = false;
            end
            
            obj = obj@DirectTrajectoryOptimization(plant,N,duration,options);
        end
        
        function obj = setupVariables(obj, N)
            switch obj.options.integration_method
                case VariationalTrajectoryOptimization.MIDPOINT
                    
                    obj.integration_method = VariationalTrajectoryOptimization.MIDPOINT;
                    
                    [phi,~,d] = obj.plant.contactConstraints(getZeroConfiguration(obj.plant), obj.options.multiple_contacts, obj.options.active_collision_options);
                    
                    %Find indices that correspond to joint angles so we can
                    %watch out for wrap-around issues
                    obj.angle_inds = [obj.plant.body(~[obj.plant.body.floating]).position_num]';
                    obj.angle_inds = obj.angle_inds(obj.angle_inds>0);
                    
                    nH = N-1;
                    nQ = obj.plant.getNumPositions();
                    nU = obj.plant.getNumInputs();
                    nC = length(phi);
                    nD = 2*length(d);
                    nL = 1+(2+2*nD)*nC; %Includes contact forces, friction cone, tangential speed, and smoothing parameter
                    
                    obj.N = N;
                    obj.nC = nC;
                    obj.nD = nD;
                    obj.nL = nL;
                    
                    num_vars = nH + nQ + N*nQ + (N-1)*nU + (N-1)*nL;
                    
                    obj.h_inds = (1:nH)';
                    obj.v0_inds = nH + (1:nQ)';
                    obj.x_inds = reshape(nH + nQ + (1:(N*nQ)), nQ, N);
                    obj.u_inds = reshape(nH + nQ + N*nQ + (1:((N-1)*nU)), nU, N-1);
                    obj.psi_inds = reshape(nH + nQ + N*nQ+(N-1)*nU + (1:((N-1)*nC)), nC, N-1);
                    obj.eta_inds = reshape(nH + nQ + N*nQ + (N-1)*nU + (N-1)*nC + (1:((N-1)*nD*nC)), nD*nC, N-1);
                    obj.c_inds = reshape(nH + nQ + N*nQ + (N-1)*nU + (N-1)*nC + (N-1)*nD*nC + (1:((N-1)*nC)), nC, N-1);
                    obj.b_inds = reshape(nH + nQ + N*nQ + (N-1)*nU + 2*(N-1)*nC + (N-1)*nD*nC + (1:((N-1)*nC*nD)), nD*nC, N-1);
                    obj.s_inds = nH + nQ + N*nQ + (N-1)*nU + 2*(N-1)*nC + 2*(N-1)*nD*nC + (1:(N-1))';
                    
                    x_names = cell(num_vars,1);
                    for j = 1:nQ
                        x_names{nH + j}=sprintf('v0[%d]',j);
                    end
                    for i = 1:N
                        for j = 1:nQ
                            x_names{nH + nQ + (i-1)*nQ+j}=sprintf('q%d[%d]',j,i);
                        end
                        if(i<N)
                            x_names{i} = sprintf('h[%d]',i);
                            for j = 1:nU
                                x_names{nH + nQ + N*nQ + (i-1)*nU+j} = sprintf('u%d[%d]',j,i);
                            end
                            for j = 1:nC
                                x_names{nH + nQ + N*nQ + (N-1)*nU + (i-1)*nC+j} = sprintf('psi%d[%d]',j,i);
                                x_names{nH + nQ + N*nQ + (N-1)*nU + (N-1)*nC + (N-1)*nC*nD + (i-1)*nC+j} = sprintf('c%d[%d]',j,i);
                            end
                            for j = 1:(nC*nD)
                                x_names{nH + nQ + N*nQ + (N-1)*nU + (N-1)*nC + (i-1)*nC*nD+j} = sprintf('eta%d[%d]',j,i);
                                x_names{nH + nQ + N*nQ + (N-1)*nU + 2*(N-1)*nC + (N-1)*nC*nD + (i-1)*nC*nD+j} = sprintf('b%d[%d]',j,i);
                            end
                            x_names{nH + nQ + N*nQ + (N-1)*nU + 2*(N-1)*nC + 2*(N-1)*nC*nD + i} = sprintf('s[%d]',i);
                        end
                    end
                    
                    obj = obj.addDecisionVariable(num_vars,x_names);
                    
                case VariationalTrajectoryOptimization.SIMPSON
                    obj.IntegrationMethod = VariationalTrajectoryOptimization.SIMPSON;
                    
                    
            end
        end
        
        function obj = addDynamicConstraints(obj)
            nQ = obj.plant.getNumPositions();
            nV = obj.plant.getNumVelocities();
            nU = obj.plant.getNumInputs();
            nC = obj.nC;
            nD = obj.nD;
            N = obj.N;
            assert(nQ == nV) % can't handle quaternions yet
            
            switch obj.integration_method
                case VariationalTrajectoryOptimization.MIDPOINT
                    
                    cnstr_opts.grad_level = 1;
                    cnstr_opts.grad_method = 'user';
                    %                     cnstr_opts.grad_method = 'numerical';
                    %                     cnstr_opts.diff_type = 'central';
                    
                    dyn_constraints = cell(N-2,1);
                    dyn_inds = cell(N-2,1);
                    
                    cont_constraints = cell(N-1,1);
                    cont_inds = cell(N-1,1);
                    
                    ineq_constraints = cell(N-1,1);
                    ineq_inds = cell(N-1,1);
                    
                    s_constraints = cell(N-1,1);
                    
                    nvars1 = 2+3*nQ+2*nU+nC+nD*nC;
                    cnstr1 = FunctionHandleConstraint(zeros(nQ,1), zeros(nQ,1), nvars1, @obj.midpoint_dynamics_constraint_fun, cnstr_opts);
                    
                    if nC
                        nvars2 = 1+2*nQ+2*nC+2*nD*nC+1;
                        %nD*nC =, 2*nC >=, 3 <=
                        cnstr2 = FunctionHandleConstraint([zeros(nD*nC+2*nC,1); -inf*ones(3,1)], [zeros(nD*nC,1); inf*ones(2*nC,1); zeros(3,1)], nvars2, @obj.midpoint_contact_constraint_fun, cnstr_opts);
                        cnstr3 = BoundingBoxConstraint(zeros(obj.nL-1,1), inf*ones(obj.nL-1,1));
                        cnstr4 = BoundingBoxConstraint(1e-6, Inf);
                        
                        for i = 1:obj.N-1
                            cont_inds{i} = {obj.h_inds(i); obj.x_inds(:,i); obj.x_inds(:,i+1); obj.psi_inds(:,i); obj.eta_inds(:,i); obj.c_inds(:,i); obj.b_inds(:,i); obj.s_inds(i)};
                            cont_constraints{i} = cnstr2.setName(sprintf('contact[%d]',i));
                            obj = obj.addConstraint(cont_constraints{i}, cont_inds{i});
                            
                            ineq_inds{i} = {[obj.psi_inds(:,i); obj.eta_inds(:,i); obj.c_inds(:,i); obj.b_inds(:,i)]};
                            ineq_constraints{i} = cnstr3.setName(sprintf('positivity[%d]',i));
                            obj = obj.addConstraint(ineq_constraints{i}, ineq_inds{i});
                            
                            s_constraints{i} = cnstr4.setName(sprintf('s[%d]',i));
                            obj = obj.addConstraint(s_constraints{i}, obj.s_inds(i));
                        end
                        obj = obj.addCost(FunctionHandleObjective(length(obj.s_inds),@(s)scost(obj,s),1), obj.s_inds(:));
                        
                        if obj.options.add_ccost
                            for ii=1:obj.N-2
                                obj = obj.addCost(FunctionHandleObjective(obj.nC*2, @(c1,c2)ccost(obj,c1,c2)),{obj.c_inds(:,ii);obj.c_inds(:,ii+1)});
                            end
                        end
                    end
                    
                    dyn_inds{1} = {obj.h_inds(1); obj.x_inds(:,1); obj.v0_inds(:); obj.x_inds(:,2); obj.u_inds(:,1); obj.c_inds(:,1); obj.b_inds(:,1)};
                    dyn_constraints{1} = FunctionHandleConstraint(zeros(nQ,1), zeros(nQ,1), 1+3*nQ+nU+nC+nD*nC, @obj.midpoint_first_step_fun, cnstr_opts);
                    dyn_constraints{1} = dyn_constraints{1}.setName('dynamics[0]');
                    obj = obj.addConstraint(dyn_constraints{1}, dyn_inds{1});
                    
                    for i = 2:obj.N-1
                        dyn_inds{i} = {obj.h_inds(i-1); obj.h_inds(i); obj.x_inds(:,i-1); obj.x_inds(:,i); obj.x_inds(:,i+1); obj.u_inds(:,i-1); obj.u_inds(:,i); obj.c_inds(:,i); obj.b_inds(:,i)};
                        dyn_constraints{i} = cnstr1.setName(sprintf('dynamics[%d]',i));
                        obj = obj.addConstraint(dyn_constraints{i}, dyn_inds{i});
                    end
                    
                case VariationalTrajectoryOptimization.SIMPSON
                    
            end
            
        end
        
        function obj = addStateConstraint(obj,constraint,time_index,x_indices)
            error('Not implemented yet');
        end
        
        function obj = addPositionConstraint(obj,constraint,time_index,x_indices)
            if ~iscell(time_index)
                % then use { time_index(1), time_index(2), ... } ,
                % aka independent constraints for each time
                time_index = num2cell(reshape(time_index,1,[]));
            end
            if nargin < 4
                x_indices = 1:size(obj.x_inds,1);
            end
            
            for j=1:length(time_index)
                
                cstr_inds = mat2cell(obj.x_inds(x_indices,time_index{j}),numel(x_indices),ones(1,length(time_index{j})));
                
                % record constraint for posterity
                obj.constraints{end+1}.constraint = constraint;
                obj.constraints{end}.var_inds = cstr_inds;
                obj.constraints{end}.time_index = time_index;
                
                obj = obj.addConstraint(constraint,cstr_inds);
            end
        end
        
        function obj = addVelocityConstraint(obj,constraint,time_index,x_indices)
            if ~iscell(time_index)
                % then use { time_index(1), time_index(2), ... } ,
                % aka independent constraints for each time
                time_index = num2cell(reshape(time_index,1,[]));
            end
            if nargin < 4
                x_indices = 1:size(obj.x_inds,1);
            end
            
            for j=1:length(time_index)
                
                if time_index{j} == 1
                    cstr_inds = mat2cell(obj.v0_inds(x_indices),numel(x_indices),ones(1,length(time_index{j})));
                else
                    error('Not implemented yet');
                    %cstr_inds = mat2cell(obj.x_inds(x_indices,time_index{j}),numel(x_indices),ones(1,length(time_index{j})));
                end
                
                % record constraint for posterity
                obj.constraints{end+1}.constraint = constraint;
                obj.constraints{end}.var_inds = cstr_inds;
                obj.constraints{end}.time_index = time_index;
                
                obj = obj.addConstraint(constraint,cstr_inds);
            end
        end
        
        function obj = addRunningCost(obj,running_cost_function)
            nQ = obj.plant.getNumPositions();
            nU = obj.plant.getNumInputs();
            switch obj.integration_method
                case VariationalTrajectoryOptimization.MIDPOINT
                    for i = 1:obj.N-1
                        running_cost = FunctionHandleObjective(1+2*nQ+nU, @(h,q1,q2,u)midpoint_running_fun(obj,running_cost_function,h,q1,q2,u));
                        inds_i = {obj.h_inds(i); obj.x_inds(:,i); obj.x_inds(:,i+1); obj.u_inds(:,i)};
                        obj = obj.addCost(running_cost,inds_i);
                    end
                case VariationalTrajectoryOptimization.SIMPSON
                    
            end
        end
        
        function obj = addNormalForceCost(obj,cost_function)
            nC = obj.nC;
            switch obj.integration_method
                case VariationalTrajectoryOptimization.MIDPOINT
                    for i = 1:obj.N-1
                        contact_cost = FunctionHandleObjective(nC, @(c)cost_function(c));
                        inds = {obj.c_inds(:,i)};
                        obj = obj.addCost(contact_cost,inds);
                    end
                case VariationalTrajectoryOptimization.SIMPSON
            end
        end
        
        function obj = addInitialCost(obj,initial_cost)
            
        end
        
        function obj = addFinalCost(obj,final_cost_function)
            nQ = obj.plant.getNumPositions();
            switch obj.integration_method
                case VariationalTrajectoryOptimization.MIDPOINT
                    final_cost = FunctionHandleObjective(2+2*nQ, @(h,q1,q2)midpoint_final_fun(obj,final_cost_function,h,q1,q2));
                    inds_i = {obj.h_inds(:); obj.x_inds(:,end-1); obj.x_inds(:,end)};
                    obj = obj.addCost(final_cost,inds_i);
                case VariationalTrajectoryOptimization.SIMPSON
                    
            end
        end
        
        function [f,df] = midpoint_first_step_fun(obj,h,q0,v0,q1,u,c,b)
            
            xin = [h;q0;v0;q1;u;c;b];
            [f,df] = midpoint_first_step(obj,xin);
            
            %             df_fd = zeros(size(df));
            %             dxin = 1e-6*eye(length(xin));
            %             for k = 1:length(xin)
            %                 df_fd(:,k) = (midpoint_first_step(obj,xin+dxin(:,k)) - midpoint_first_step(obj,xin-dxin(:,k)))/2e-6;
            %             end
            %
            %             disp('First step derivative error:');
            %             disp(max(abs(df_fd(:)-df(:))));
            
        end
        
        
        function [f,df] = midpoint_first_step(obj,xin)
            
            nC = obj.nC;
            nD = obj.nD;
            nQ = obj.plant.getNumPositions();
            nU = obj.plant.getNumInputs();
            
            h = xin(1);
            q0 = xin(1+(1:nQ));
            v0 = xin(1+nQ+(1:nQ));
            q1 = xin(1+2*nQ+(1:nQ));
            u = xin(1+3*nQ+(1:nU));
            c = xin(1+3*nQ+nU+(1:nC));
            b = xin(1+3*nQ+nU+nC+(1:nD*nC));
            
            %Discrete Euler-Lagrange equation
            [M0,~,~,dM0] = manipulatorDynamics(obj.plant, q0, zeros(nQ,1));
            dM0 = reshape(dM0,nQ*nQ,2*nQ);
            dM0 = dM0(:,1:nQ);
            p0 = M0*v0;
            qm = obj.qavg(q0,q1);
            vm = obj.qdiff(q0,q1,h);
            [D1L,D2L,D1D1L,D1D2L,D2D2L,B,dB] = obj.LagrangianDerivs(qm,vm);
            f_del = p0 + (h/2)*D1L - D2L;
            
            df_del = [(1/2)*D1L - ((h/2)*D1D2L'-D2D2L)*vm/h, ... % d/dh
                kron(v0',eye(nQ))*dM0 + (h/4)*D1D1L - (1/2)*D1D2L' - (1/2)*D1D2L + (1/h)*D2D2L, ... % d/dq0
                M0, ... % d/dv0
                (h/4)*D1D1L + (1/2)*D1D2L' - (1/2)*D1D2L - (1/h)*D2D2L, ... % d/dq1
                zeros(nQ,nU), zeros(nQ,nC), zeros(nQ,nD*nC)]; % d/du, d/dc, d/db
            
            %Contact basis
            kinopts = struct();
            kinopts.compute_gradients = true;
            kin = obj.plant.doKinematics(q1, vm, kinopts);
            [~,~,~,~,~,~,~,~,n,D,dn,dD] = obj.plant.contactConstraints(kin, obj.options.multiple_contacts,obj.options.active_collision_options);
            if isempty(n)
                n = zeros(0,nQ);
                dn = zeros(0,nQ);
            end
            D = reshape(cell2mat(D')',nQ,nC*nD)';
            dD = reshape(cell2mat(dD)',nQ,nC*nD*nQ)';
            
            f = f_del + (h/2)*B*u + h*(n'*c + D'*b);
            
            df = df_del + [(1/2)*B*u + n'*c + D'*b, ... % d/dh
                (h/4)*kron(u',eye(nQ))*dB, ... % d/dq0
                zeros(nQ,nQ), ... % d/dv0
                (h/4)*kron(u',eye(nQ))*dB + h*kron(c',eye(nQ))*comm(nC,nQ)*dn + h*kron(b',eye(nQ))*comm(nD*nC,nQ)*dD, ... % d/dq1
                (h/2)*B, h*n', h*D']; % d/du, d/dc, d/db
        end
        
        function [f,df] = midpoint_dynamics_constraint_fun(obj,h1,h2,q1,q2,q3,u1,u2,c2,b2)
            
            xin = [h1;h2;q1;q2;q3;u1;u2;c2;b2];
            [f,df] = obj.midpoint_dynamics(xin);
%             fprintf('Dynamics: %f \r', max(abs(f)));
            
            %             df_fd = zeros(size(df));
            %             step = sqrt(eps(max(xin)));
            %             dxin = step*eye(length(xin));
            %             for k = 1:length(xin)
            %                 df_fd(:,k) = obj.qdiff(obj.midpoint_dynamics(xin-dxin(:,k)), obj.midpoint_dynamics(xin+dxin(:,k)), 2*step);
            %             end
            
            %             disp('Dynamics Derivative Error:');
            %             disp(max(abs(df_fd(:)-df(:))));
            
        end
        
        function [f,df] = midpoint_dynamics(obj,xin)
            
            nC = obj.nC;
            nD = obj.nD;
            nQ = obj.plant.getNumPositions();
            nU = obj.plant.getNumInputs();
            
            h1 = xin(1);
            h2 = xin(2);
            q1 = xin(2+(1:nQ));
            q2 = xin(2+nQ+(1:nQ));
            q3 = xin(2+2*nQ+(1:nQ));
            u1 = xin(2+3*nQ+(1:nU));
            u2 = xin(2+3*nQ+nU+(1:nU));
            c2 = xin(2+3*nQ+2*nU+(1:nC));
            b2 = xin(2+3*nQ+2*nU+nC+(1:nC*nD));
            
            %Take care of angle wrap-around
            qm1 = obj.qavg(q1,q2);
            vm1 = obj.qdiff(q1,q2,h1);
            qm2 = obj.qavg(q2,q3);
            vm2 = obj.qdiff(q2,q3,h2);
            
            %Discrete Euler-Lagrange equation
            [D1L1,D2L1,D1D1L1,D1D2L1,D2D2L1,B1,dB1] = obj.LagrangianDerivs(qm1,vm1);
            [D1L2,D2L2,D1D1L2,D1D2L2,D2D2L2,B2,dB2] = obj.LagrangianDerivs(qm2,vm2);
            f_del = (h1/2)*D1L1 + D2L1 + (h2/2)*D1L2 - D2L2;
            
            df_del = [0.5*D1L1 - ((h1/2)*D1D2L1'+D2D2L1)*(vm1/h1), 0.5*D1L2 - ((h2/2)*D1D2L2'-D2D2L2)*(vm2/h2), ... % d/dh1, d/dh2
                (h1/2)*((1/2)*D1D1L1-(1/h1)*D1D2L1')+(1/2)*D1D2L1-(1/h1)*D2D2L1, ... % d/dq1
                (h1/2)*((1/2)*D1D1L1+(1/h1)*D1D2L1')+(1/2)*D1D2L1+(1/h1)*D2D2L1 + (h2/2)*((1/2)*D1D1L2-(1/h2)*D1D2L2')-(1/2)*D1D2L2+(1/h2)*D2D2L2, ... % d/dq2
                (h2/2)*((1/2)*D1D1L2+(1/h2)*D1D2L2') - (1/2)*D1D2L2 - (1/h2)*D2D2L2, ... % d/dq3
                zeros(nQ, 2*nU+nC+(nC*nD))];
            
            %Contact basis
            kinopts = struct();
            kinopts.compute_gradients = true;
            kin2 = obj.plant.doKinematics(q3, vm2, kinopts);
            [~,~,~,~,~,~,~,~,n2,D2,dn2,dD2] = obj.plant.contactConstraints(kin2, obj.options.multiple_contacts,obj.options.active_collision_options);
            if isempty(n2)
                n2 = zeros(0,nQ);
                dn2 = zeros(0,nQ);
            end
            D2 = reshape(cell2mat(D2')',nQ,nC*nD)';
            dD2 = reshape(cell2mat(dD2)',nQ,nC*nD*nQ)';
            
            %Total dynamics residual incluing control + contact forces
            f = f_del + (h1/2)*B1*u1 + (h2/2)*B2*u2 + h2*(n2'*c2 + D2'*b2);
            
            %Dynamics Derivatives
            df = df_del + [(1/2)*B1*u1, (1/2)*B2*u2 + n2'*c2 + D2'*b2, ... % d/dh1, d/dh2
                (h1/4)*kron(u1', eye(nQ))*dB1, ... % d/dq1
                (h1/4)*kron(u1', eye(nQ))*dB1 + (h2/4)*kron(u2', eye(nQ))*dB2, ... % d/dq2
                (h2/4)*kron(u2', eye(nQ))*dB2 + h2*kron(c2', eye(nQ))*comm(nC,nQ)*dn2 + h2*kron(b2', eye(nQ))*comm(nD*nC,nQ)*dD2, ... % d/dq3
                (h1/2)*B1, (h2/2)*B2, h2*n2', h2*D2']; % d/du1, d/du2, d/dc2, d/db2
        end
        
        function [f,df] = midpoint_contact_constraint_fun(obj,h,q1,q2,psi,eta,c,b,s)
            
            xin = [h;q1;q2;psi;eta;c;b;s];
            [f,df] = obj.midpoint_contact(xin);
%             fprintf('Slack: %f \r', s);
            
            %             df_fd = zeros(size(df));
            %             dxin = 1e-6*eye(length(xin));
            %             for k = 1:length(xin)
            %                 df_fd(:,k) = (obj.midpoint_contact(xin+dxin(:,k)) - obj.midpoint_contact(xin-dxin(:,k)))/2e-6;
            %             end
            %
            %             disp('Contact Derivative Error:');
            %             disp(max(abs(df_fd(:)-df(:))));
            
        end
        
        function [f,df] = midpoint_contact(obj,xin)
            mu = 1; %This is currently hard-coded in Drake
            nC = obj.nC;
            nD = obj.nD;
            nQ = obj.plant.getNumPositions();
            
            h = xin(1);
            q1 = xin(1+(1:nQ));
            q2 = xin(1+nQ+(1:nQ));
            psi = xin(1+2*nQ+(1:nC));
            eta = xin(1+2*nQ+nC+(1:nD*nC));
            c = xin(1+2*nQ+nC+nD*nC+(1:nC));
            b = xin(1+2*nQ+2*nC+nD*nC+(1:nD*nC));
            s = xin(end);
            
            vm = obj.qdiff(q1,q2,h);
            
            %Contact basis
            kinopts = struct();
            kinopts.compute_gradients = true;
            kin = obj.plant.doKinematics(q2, vm, kinopts);
            [phi,~,~,~,~,~,~,~,n,D,dn,dD] = obj.plant.contactConstraints(kin, obj.options.multiple_contacts, obj.options.active_collision_options);
            if isempty(n)
                n = zeros(0,nQ);
                dn = zeros(0,nQ);
            end
            D = reshape(cell2mat(D')',nQ,nC*nD)';
            dD = reshape(cell2mat(dD)',nQ,nC*nD*nQ)';
            
            e = ones(nD,1);
            E = kron(eye(nC),e');
            
            %Tangential velocity
            f1 = D*vm + E'*psi - eta; % = 0
            
            %Signed distance
            g1 = phi; % >= 0
            
            %Friction cone
            g2 = mu*c - E*b; % >= 0
            
            %Normal force complementarity
            l1 = phi'*c - s; % <= 0
            
            %Tangential velocity complementarity
            l2 = (mu*c - E*b)'*psi - s; % <= 0
            
            %Friction complementarity
            l3 = eta'*b - s; % <= 0
            
            f = [f1; g1; g2; l1; l2; l3];
            
            %xin = [h;q1;q2;psi;eta;c;b;s];
            df = [-D*vm/h, -D/h, D/h + kron(vm', eye(nD*nC))*dD, E', -eye(nD*nC), zeros(nD*nC,nC), zeros(nD*nC,nD*nC), zeros(nD*nC,1);
                zeros(nC,1), zeros(nC,nQ), n, zeros(nC,nC), zeros(nC,nD*nC), zeros(nC,nC), zeros(nC,nD*nC), zeros(nC,1);
                zeros(nC,1), zeros(nC,nQ), zeros(nC,nQ), zeros(nC,nD*nC), zeros(nC,nC), mu*eye(nC), -E, zeros(nC,1);
                zeros(1,1), zeros(1,nQ), c'*n, zeros(1,nC), zeros(1,nD*nC), phi', zeros(1,nD*nC), -1;
                zeros(1,1), zeros(1,nQ), zeros(1,nQ), (mu*c - E*b)', zeros(1,nD*nC), psi'*mu, -psi'*E, -1;
                zeros(1,1), zeros(1,nQ), zeros(1,nQ), zeros(1,nC), b', zeros(1,nC), eta', -1];
        end
        
        function [D1L,D2L,D1D1L,D1D2L,D2D2L,B,dBdq] = LagrangianDerivs(obj,q,v)
            nq = length(q);
            nv = length(v);
            [M,G,B,dM,dG,dB] = manipulatorDynamics(obj.plant, q, zeros(nv,1));
            
            dM = reshape(dM,nq*nq,nq+nv);
            dMdq = dM(:,1:nq);
            dBdq = dB(:,1:nq);
            
            D1L = 0.5*dMdq'*kron(v,v) - G;
            D2L = M*v;
            
            %D1D1L = -dG(:,1:nq); %throwing out second derivative of M terms here
            
            D1D1L = zeros(nq);
            step = sqrt(eps(max(q)));
            deltaq = step*eye(nq);
            for k = 1:nq
                [~,Gp,~,dMp] = manipulatorDynamics(obj.plant, q+deltaq(:,k), zeros(nv,1));
                dMp = reshape(dMp,nq*nq,nq+nv);
                dMdqp = dMp(:,1:nq);
                
                [~,Gm,~,dMm] = manipulatorDynamics(obj.plant, q-deltaq(:,k), zeros(nv,1));
                dMm = reshape(dMm,nq*nq,nq+nv);
                dMdqm = dMm(:,1:nq);
                
                D1p = 0.5*dMdqp'*kron(v,v) - Gp;
                D1m = 0.5*dMdqm'*kron(v,v) - Gm;
                
                D1D1L(:,k) = (D1p - D1m)/(2*step);
            end
            
            %disp(sprintf('D1D1L error: %d',max(abs(D1D1L_fd(:)-D1D1L(:)))));
            
            D1D2L = kron(v',eye(nq))*dMdq;
            D2D2L = M;
            
        end
        
        function [c,dc] = ccost(~,c1,c2)
            cdiff = c1-c2;
            c = 0.5*(cdiff'*cdiff);
            I = eye(length(c1));
            dc = [cdiff'*I,-cdiff'*I];
        end
        
        function [c,dc] = scost(obj, s)
            c = obj.options.s_weight*sum(s);
            dc = obj.options.s_weight*ones(1,length(s));
        end
        
        function [xtraj,utraj,ctraj,btraj,psitraj,etatraj,straj,z,F,info,infeasible_constraint_name] = solveTraj(obj,t_init,traj_init)
            [xtraj,utraj,z,F,info,infeasible_constraint_name] = solveTraj@DirectTrajectoryOptimization(obj,t_init,traj_init);
            t = [0; cumsum(z(obj.h_inds))];
            if obj.nC>0
                ctraj = PPTrajectory(zoh(t,[reshape(z(obj.c_inds),[],obj.N-1),z(obj.c_inds(:,end))]));
                btraj = PPTrajectory(zoh(t,[reshape(z(obj.b_inds),[],obj.N-1),z(obj.b_inds(:,end))]));
                psitraj = PPTrajectory(zoh(t,[reshape(z(obj.psi_inds),[],obj.N-1),z(obj.psi_inds(:,end))]));
                etatraj = PPTrajectory(zoh(t,[reshape(z(obj.eta_inds),[],obj.N-1),z(obj.eta_inds(:,end))]));
                straj = PPTrajectory(zoh(t,[reshape(z(obj.s_inds),[],obj.N-1),z(obj.s_inds(end))]));
            else
                ctraj = [];
                btraj = [];
                psitraj = [];
                etatraj = [];
                straj = [];
            end
        end
        
        function [f,df] = midpoint_running_fun(obj,running_fun,h,q1,q2,u)
            xin = [h;q1;q2;u];
            [f,df] = midpoint_running(obj,running_fun,xin);
            
            %             df_fd = zeros(size(df));
            %             deltax = 1e-6*eye(length(xin));
            %             for k = 1:length(xin)
            %                 df_fd(:,k) = (midpoint_running(obj,running_fun,xin+deltax(:,k)) - midpoint_running(obj,running_fun,xin-deltax(:,k)))/2e-6;
            %             end
            %
            %             disp('Running cost derivative error:');
            %             disp(max(abs(df_fd(:)-df(:))));
        end
        
        function [f,df] = midpoint_running(obj,running_fun,xin)
            
            nQ = obj.plant.getNumPositions();
            
            h = xin(1);
            q1 = xin(1+(1:nQ));
            q2 = xin(1+nQ+(1:nQ));
            u = xin((2+2*nQ):end);
            
            qm = obj.qavg(q1,q2);
            vm = obj.qdiff(q1,q2,h);
            
            [f,dg] = running_fun(h,[qm; vm],u);
            df = [dg(:,1)-dg(:,(1+nQ)+(1:nQ))*vm/h, 0.5*dg(:,1+(1:nQ))-(1/h)*dg(:,(1+nQ)+(1:nQ)), 0.5*dg(:,1+(1:nQ))+(1/h)*dg(:,(1+nQ)+(1:nQ)), dg(:,(2+2*nQ):end)];
        end
        
        function [f,df] = midpoint_final_fun(obj,final_fun,h,q1,q2)
            nQ = obj.plant.getNumPositions();
            tf = sum(h);
            vm = obj.qdiff(q1,q2,h(end));
            [f,dg] = final_fun(tf,[q2; vm]);
            df = [kron(ones(1,obj.N-2),dg(:,1)), dg(:,1)-(dg(:,(1+nQ)+(1:nQ))*vm/h(end)), -(1/h(end))*dg(:,nQ+(1:nQ)), dg(:,1:nQ)+(1/h(end))*dg(:,nQ+(1:nQ))];
        end
        
        function xtraj = reconstructStateTrajectory(obj,z)
            nQ = obj.plant.getNumPositions();
            nC = obj.nC;
            nD = obj.nD;
            
            t = [0; cumsum(z(obj.h_inds))];
            x = reshape(z(obj.x_inds),[],obj.N);
            
            switch obj.integration_method
                case VariationalTrajectoryOptimization.MIDPOINT
                    
                    q = x(1:nQ,:);
                    qtraj = PPTrajectory(foh(t,q));
                    
                    v = [diff(x,1,2)/z(obj.h_inds(1)), zeros(nQ,1)]; %zoh (correctly in this case) ignores the last entry in v
                    vtraj = PPTrajectory(zoh(t,v));
                    
                    %                     h = z(obj.h_inds);
                    %                     u = z(obj.u_inds);
                    %                     c = z(obj.c_inds);
                    %                     b = z(obj.b_inds);
                    %
                    %                     v = zeros(size(q));
                    %                     v(:,1) = z(obj.v0_inds);
                    %                     for k = 1:(obj.N-1)
                    %                         qm = obj.qavg(q(:,k),q(:,k+1));
                    %                         vm = obj.qdiff(q(:,k),q(:,k+1),h(k));
                    %
                    %                         [D1L,D2L,~,~,~,B] = obj.LagrangianDerivs(qm,vm);
                    %                         M1 = manipulatorDynamics(obj.plant, q(:,k+1), vm);
                    %
                    %                         %Contact basis
                    %                         kin = obj.plant.doKinematics(q(:,k+1));
                    %                         [~,~,~,~,~,~,~,~,n,D] = obj.plant.contactConstraints(kin, obj.options.multiple_contacts,obj.options.active_collision_options);
                    %                         if isempty(n)
                    %                             n = zeros(0,nQ);
                    %                         end
                    %                         D = reshape(cell2mat(D')',nQ,nC*nD)';
                    %
                    %                         p1 = (h(k)/2)*D1L + D2L + (h(k)/2)*B*u(:,k) + (h(k)/2)*(n'*c(:,k) + D'*b(:,k));
                    %
                    %                         v(:,k+1) = M1\p1;
                    %                     end
                    %
                    %                     vtraj = PPTrajectory(foh(t,v));
                    
                case VariationalTrajectoryOptimization.SIMPSON
                    
            end
            
            xtraj = [qtraj; vtraj];
            xtraj = xtraj.setOutputFrame(obj.plant.getStateFrame);
        end
        
        function utraj = reconstructInputTrajectory(obj,z)
            if size(obj.u_inds,1) > 0
                nU = obj.plant.getNumInputs();
                switch obj.integration_method
                    case VariationalTrajectoryOptimization.MIDPOINT
                        t = [0; cumsum(z(obj.h_inds))];
                        u = [reshape(z(obj.u_inds),[],obj.N-1), zeros(nU,1)]; %zoh (correctly in this case) ignores the last entry in v
                        utraj = PPTrajectory(zoh(t,u));
                        utraj = utraj.setOutputFrame(obj.plant.getInputFrame);
                    case VariationalTrajectoryOptimization.SIMPSON
                        
                end
            else
                utraj=[];
            end
        end
        
        function qm = qavg(obj,q1,q2)
            qm = (q1+q2)/2;
            if ~isempty(obj.angle_inds)
                qm(obj.angle_inds) = angleAverage(q1(obj.angle_inds),q2(obj.angle_inds));
            end
        end
        
        function vm = qdiff(obj,q1,q2,h)
            vm = (q2-q1)/h;
            vm(obj.angle_inds) = angleDiff(q1(obj.angle_inds),q2(obj.angle_inds))/h;
        end
        
        function z0 = getInitialVars(obj,t_init,traj_init)
            nQ = obj.plant.getNumPositions();
            
            if isscalar(t_init)
                t_init = linspace(0,t_init,obj.N);
            elseif length(t_init) ~= obj.N
                error('The initial sample times must have the same length as property N')
            end
            z0 = zeros(obj.num_vars,1);
            z0(obj.h_inds) = diff(t_init);
            
            if nargin < 3
                traj_init = struct();
            end
            
            nU = getNumInputs(obj.plant);
            if isfield(traj_init,'u')
                z0(obj.u_inds) = traj_init.u.eval(t_init(1:end-1));
            else
                z0(obj.u_inds) = 0.01*randn(nU,obj.N-1);
            end
            
            if isfield(traj_init,'x')
                xsamp = traj_init.x.eval(t_init);
                z0(obj.x_inds) = xsamp(1:nQ,:);
                z0(obj.v0_inds) = xsamp(nQ+(1:nQ),1);
            end
            
            if isfield(traj_init,'c')
                z0(obj.c_inds) = traj_init.c.eval(t_init(1:end-1));
            else
                z0(obj.c_inds(:)) = .1*rand(length(obj.c_inds(:)),1);
            end
            
            if isfield(traj_init,'b')
                z0(obj.b_inds) = traj_init.b.eval(t_init(1:end-1));
            else
                z0(obj.b_inds(:)) = .1*rand(length(obj.b_inds(:)),1);
            end
            
            if isfield(traj_init,'psi')
                z0(obj.psi_inds) = traj_init.psi.eval(t_init(1:end-1));
            else
                z0(obj.psi_inds(:)) = .1*rand(length(obj.psi_inds(:)),1);
            end
            
            if isfield(traj_init,'eta')
                z0(obj.eta_inds) = traj_init.eta.eval(t_init(1:end-1));
            else
                z0(obj.eta_inds(:)) = .1*rand(length(obj.eta_inds(:)),1);
            end
            
            if isfield(traj_init,'s')
                z0(obj.s_inds) = traj_init.s.eval(t_init(1:end-1));
            else
                z0(obj.s_inds(:)) = 1;
            end
            %              temp = load('z0.mat');
            %              z0 = temp.z0;
            %              z0(obj.u_inds) = z0(obj.u_inds) + 0.01*randn(size(obj.u_inds));
            %              z0(obj.c_inds) = z0(obj.c_inds) + 0.01*randn(size(obj.c_inds));
            %              z0(obj.b_inds) = z0(obj.b_inds) + 0.01*randn(size(obj.b_inds));
            %              z0(obj.eta_inds) = z0(obj.eta_inds) + 0.01*randn(size(obj.eta_inds));
            %              z0(obj.psi_inds) = z0(obj.psi_inds) + 0.01*randn(size(obj.psi_inds));
            % %              obj.z0 = temp.z0;
        end
    end
end