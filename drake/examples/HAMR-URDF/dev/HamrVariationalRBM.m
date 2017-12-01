classdef HamrVariationalRBM < RigidBodyManipulator
    
    properties (SetAccess = protected, GetAccess = public)
        x0
        %         grav = [0; 0; -9.81e-3];
        LOOP_BODIES = {'FLS3', 'FLS_Loop', [0; 0; 0], [0.53; -0.12; 0.0];
            'FLL3', 'FLL_Loop', [0; 0; 0], [0.0; -0.12; -0.58];
            'FLS4', 'FLSFB_Loop', [0; 0; 0], [2.05236472111; -5.98079552649; 0.0];
            'RLS3', 'RLS_Loop', [0; 0; 0],[-0.53; -0.12; 0.0];
            'RLL3', 'RLL_Loop', [0; 0; 0],[0.0; -0.12; -0.58];
            'RLS4', 'RLSFB_Loop', [0; 0; 0], [-2.05236472111; -5.98079552649; 0.0];
            'FRS3', 'FRS_Loop', [0; 0; 0],[0.53; 0.12; 0.0];
            'FRL3', 'FRL_Loop', [0; 0; 0],[0.0; 0.12; -0.58];
            'FRS4', 'FRSFB_Loop', [0; 0; 0], [2.05236472111; 5.98079552649; 0.0];
            'RRS3', 'RRS_Loop', [0; 0; 0],[-0.53; 0.12; 0.0];
            'RRL3', 'RRL_Loop', [0; 0; 0],[0.0; 0.12; -0.58];
            'RRS4', 'RRSFB_Loop', [0; 0; 0], [-2.05236472111; 5.98079552649; 0.0]};
        %
        %         LOOP_BODIES = {'FLS_Loop', 'Chassis', [12.754659815390; 9.79; 1.71], [-9.39733204324347E-02; -5.48845250802898E-04; -0.379657560961688];
        %             'FLL_Loop', 'Chassis', [10.504659815390; 9.79; 0.46], [-9.39733204324347E-02; -5.48845250802898E-04; -0.379657560961688];
        %             'FLSFB_Loop', 'FLL4', [0.0; 5.98079552649; 1.64999985615], [-2.4877; 1.5376; -1.509];
        %             'RLS_Loop', 'Chassis', [-12.754659815390; 9.79; 1.71],[-9.39733204324347E-02; -5.48845250802898E-04; -0.379657560961688];
        %             'RLL_Loop', 'Chassis', [-10.504659815390; 9.79; 0.46],[-9.39733204324347E-02; -5.48845250802898E-04; -0.379657560961688];
        %             'RLSFB_Loop', 'RLL4', [0.0; 5.98079552649; 1.64999985615], [2.4877; 1.5376; -1.509];
        %             'FRS_Loop', 'Chassis', [12.754659815390; -9.79; 1.71],[-9.39733204324347E-02; -5.48845250802898E-04; -0.379657560961688];
        %             'FRL_Loop', 'Chassis', [10.504659815390; -9.79; 0.46],[-9.39733204324347E-02; -5.48845250802898E-04; -0.379657560961688];
        %             'FRSFB_Loop', 'FRL4', [0.0; -5.98079552649; 1.64999985615], [-2.4877; -1.5376; -1.509];
        %             'RRS_Loop', 'Chassis', [-12.754659815390; -9.79; 1.71],[-9.39733204324347E-02; -5.48845250802898E-04; -0.379657560961688];
        %             'RRL_Loop', 'Chassis', [-10.504659815390; -9.79; 0.46],[-9.39733204324347E-02; -5.48845250802898E-04; -0.379657560961688];
        %             'RRSFB_Loop', 'RRL4', [0.0; -5.98079552649; 1.64999985615], [2.4877; -1.5376; -1.509]};
        
    end
    
    methods
        
        function obj=HamrVariationalRBM(urdf,options)
            
            typecheck(urdf,'char');
            
            if nargin < 2
                options = struct();
            end
            
            obj = obj@RigidBodyManipulator(urdf,options);
            obj.x0 = zeros(2*obj.getNumPositions(), 1);
            
            %Change gravity
            obj = obj.setGravity([0; 0; -9.81e-3]);
            
            obj = compile(obj);
            obj = obj.addConstraints();
        end
        
        function [y, Jy] = output(obj, t, x, u)
            
            nx = numel(x);
            actuated_dof = obj.getActuatedJoints();
            y = [x; x(actuated_dof); x(actuated_dof + nx/2)];
            
            if nargout > 1
                ny = numel(y); nu = numel(u);
                Jy1 = [zeros(nx, 1), eye(nx), zeros(nx, nu)];       % derivative of state
                Jy2 = zeros(ny - nx, 1 + nx + nu);                  % derivative of actuated dof q and qd
                r_ind = (1:ny-nx)';
                c_ind = 1+[actuated_dof; actuated_dof + nx/2];
                Jy2(sub2ind([ny-nx, 1 + nx + nu], r_ind, c_ind)) = 1;
                Jy = [Jy1; Jy2];
            end
        end
        
        
        function obj = compile(obj)
            obj = compile@RigidBodyManipulator(obj);
            
            %Add Ouputs
            joint_names = {obj.body.jointname}';
            actuated_dof = obj.getActuatedJoints();
            obj = obj.setNumOutputs(obj.getNumStates()+ 2*numel(actuated_dof));
            state_frame = obj.getStateFrame();
            
            if obj.getBody(2).floating
                act_jt_names = [joint_names(actuated_dof - 4); ...
                    joint_names(actuated_dof - 4)];
                %                 act_jt_names =
            else
                act_jt_names = [joint_names(actuated_dof + 1); ...
                    joint_names(actuated_dof +1 )];
            end
            
            output_frame = MultiCoordinateFrame( ...
                {state_frame.getFrameByNum(1), ...
                state_frame.getFrameByNum(2), ...
                CoordinateFrame('ActuatorDeflectionandRate', 2*obj.getNumActuatedDOF(), ...
                {}, act_jt_names)},[ones(obj.getNumPositions(),1); ...
                2*ones(obj.getNumVelocities(),1); 3*ones(2*numel(actuated_dof),1)]);
            
            obj = obj.setOutputFrame(output_frame);
        end
        
        function [f, df, df2] = kl_const_fun(obj, x)
            
            [f,df, df2] = kl_const(obj,x);
            
%             df_fd = zeros(size(df));
%             dx = 1e-6*eye(length(x));
%             for k = 1:length(x)
%                 df_fd(:,k) = (kl_const(obj,x+dx(:,k)) - kl_const(obj,x-dx(:,k)))/2e-6;
%             end
%             
%             disp('Kinematic Loop Constraint derivative error:');
%             disp(max(abs(df_fd(:)-df(:))));
            % %
        end
        
        function [f, df, df2] = kl_const(obj, q)
            
            warning('off', 'Drake:RigidBodyManipulator:WeldedLinkInd')
            nQ = obj.getNumPositions();
            nKL = size(obj.LOOP_BODIES,1);
            
            f = zeros(6*nKL, 1);
            df = zeros(6*nKL, nQ);
            
            if nargout > 2
                df2 = zeros(6*nKL, nQ^2);
                kinsol = doKinematics(obj,q,zeros(nQ,1), struct('compute_gradients', true));
                for i = 1:nKL
                    loop_body1 = obj.findLinkId(obj.LOOP_BODIES{i,1});
                    loop_body2 = obj.findLinkId(obj.LOOP_BODIES{i,2});
                    
                    [x1, J1, dJ1] = obj.forwardKin(kinsol, loop_body1, obj.LOOP_BODIES{i,3}, struct('rotation_type', 1));
                    [x2, J2, dJ2] = obj.forwardKin(kinsol, loop_body2, obj.LOOP_BODIES{i,4}, struct('rotation_type', 1));
                    %                     [x1, J1, dJ1] = obj.forwardKin(kinsol, loop_body1, -obj.LOOP_BODIES{i,3}, struct('rotation_type', 1));
                    %                     [x2, J2, dJ2] = obj.forwardKin(kinsol, loop_body2, zeros(3,1), struct('rotation_type', 1));
                    
                    f(6*(i-1)+1:6*i) = x2 - x1;
                    df(6*(i-1)+1:6*i, :) = J2 - J1;
                    df2(6*(i-1)+1:6*i, :) = dJ2 - dJ1;
                end
            else
                kinsol = doKinematics(obj,q,zeros(nQ,1), struct('compute_gradients', false));
                for i = 1:nKL
                    loop_body1 = obj.findLinkId(obj.LOOP_BODIES{i,1});
                    loop_body2 = obj.findLinkId(obj.LOOP_BODIES{i,2});
                    
                    [x1, J1] = obj.forwardKin(kinsol, loop_body1, obj.LOOP_BODIES{i,3}, struct('rotation_type', 1));
                    [x2, J2] = obj.forwardKin(kinsol, loop_body2, obj.LOOP_BODIES{i,4}, struct('rotation_type', 1));
                    %                     [x1, J1] = obj.forwardKin(kinsol, loop_body1, -obj.LOOP_BODIES{i,3}, struct('rotation_type', 1));
                    %                     [x2, J2] = obj.forwardKin(kinsol, loop_body2, zeros(3,1), struct('rotation_type', 1));
                    
                    f(6*(i-1)+1:6*i) = x2 - x1;
                    df(6*(i-1)+1:6*i, :) = J2 - J1;
                end
            end            
            warning('on', 'Drake:RigidBodyManipulator:WeldedLinkInd')           
        end
        
        
        function obj = addConstraints(obj)
            
            % Add loop constraints
            nKL = size(obj.LOOP_BODIES,1);
            nQ = obj.getNumPositions();
            
            cnstr_opts.grad_level = 2;
            cnstr_opts.grad_method = 'user';
            
            klconst = FunctionHandleConstraint(zeros(6*nKL, 1), ...
                zeros(6*nKL, 1), nQ, @(x)kl_const_fun(obj, x), cnstr_opts);
            
            const_names = cell(6*nKL, 1);
            for i = 1:nKL
                for j = 1:6
                    const_names{(i-1)*6+j} = sprintf([obj.LOOP_BODIES{i, 1}, '[%d]',], j);
                end
            end
            
            klconst = klconst.setName(const_names);
            
            % Add Constraints
            obj2 = obj.addPositionEqualityConstraint(klconst, 1:nQ);
            klconst = klconst.setEvalHandle(@(x)kl_const_fun(obj2, x));
            obj = obj.addPositionEqualityConstraint(klconst, 1:nQ);           
            
        end     
        
        function nActuatedDOF = getNumActuatedDOF(obj)
            nActuatedDOF = numel(obj.getActuatedJoints());
        end
        
        function obj = setInitialState(obj,x0)
            typecheck(x0,'double');
            sizecheck(x0,obj.getNumStates());
            obj.x0 = x0;
        end
        
        function x0 = getInitialState(obj)
%             q0 = load('q0_biased.mat');
%             q0 = q0.q0_biased;
%             x0 = [zeros(6,1); q0(1:44); zeros(6,1); q0(45:end)];
            x0 = zeros(obj.getNumStates(),1); 
%             x0(3) = 12.04;
            obj.x0 = x0;
        end
    end
    
    
    
end
