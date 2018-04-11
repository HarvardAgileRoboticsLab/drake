classdef HAMRLegTracking < DrakeSystem
    
    
    properties (SetAccess=private)
        rwact;          % cascade with actuators
        robot;          % to be controlled
        act;            % actuators
        nq;             % number of positions
        nv;             % number of velocities
        qa;             % active dof
        x_des;          % desired position/velocity
        u0;             % feed forward voltage
        CTYPE           % feed back controller type
        K;              % full state feedback matrix
    end
    
    
    methods
        function obj = HAMRLegTracking(rwact, r, act, u0, x_des, opt)
            % @param r rigid body manipulator instance
            
            input_frame = rwact.getOutputFrame();
            output_frame = rwact.getInputFrame();
            
            obj = obj@DrakeSystem(0,0,rwact.getNumOutputs,rwact.getNumInputs,true,true);
            
            obj = setInputFrame(obj,input_frame);
            obj = setOutputFrame(obj,output_frame);
            
            obj.rwact = rwact;
            obj.robot = r;
            obj.act = act; 
            obj.nq = r.getNumPositions();
            obj.nv = r.getNumVelocities();
            obj.qa = r.getActuatedJoints();
                        
            orien = zeros(numel(obj.qa), 1);
            for i = 1:numel(obj.qa)
                orien(i) = act.dummy_bender(i).orien;
            end
            
            obj.CTYPE = opt.ctype;
            
            switch obj.CTYPE
                
                case 'legpd'
                    obj.K = [diag(-orien)*opt.kp, diag(-orien)*opt.kd];
                case 'actlqr'
                    lti_approx = load('TYM_LinearEstimateDrake.mat');     
                    nx_sl = size(lti_approx.Ap, 2);
                    nu_sl = size(lti_approx.Bp, 2);
                    nl_lti = size(lti_approx.Ap,3);
                    
                    nx_lti = nx_sl*nl_lti;
                    nu_lti = nu_sl*nl_lti;                    
                    
                    Alti = zeros(nx_lti);
                    Blti = zeros(nx_lti, nu_lti);
                    
                    % unwrap A, B matrices
                    for i = 1:nl_lti
                        Alti((i-1)*nx_sl+(1:nx_sl), (i-1)*nx_sl+(1:nx_sl)) = lti_approx.Ap(:,:,i);
                        Blti((i-1)*nx_sl+(1:nx_sl), (i-1)*nu_sl+(1:nu_sl)) = lti_approx.Bp(:,:,i);
                    end
                    
%                     % LQR gains
                    e = [opt.Qpos; opt.Qvel]*ones(1, nx_lti/2);
                    Qlti = diag(e(:)); 
%                     Qlti([1,2, 5, 6, 9, 10, 13, 14], :) = 0; 
                    Rlti = opt.rho*eye(nu_lti); 
                    [P, ~, obj.K] = dare(Alti, Blti, Qlti, Rlti);

                    
                case 'act_tvlqr'
            end
            
            obj.u0 = u0;
            obj.x_des = x_des;
        end
        
        
        %         function y = output(obj, t, ~, x)
        %
        %             u0 = obj.u0.eval(t);
        %             x_des = obj.x_des.eval(t);
        %             x_hat = x([obj.qa; obj.nq+obj.qa]);
        %
        %             err = x_des - x_hat;
        %
        %             y = u0 + [obj.Kp, obj.Kd]*err;
        %
        % %             disp(
        %
        %             y(y > 0.3) = 0.3;
        %             y(y < -0.3) = -0.3;
        %
        %         end
        
        
        function u = output(obj, t, ~, x)
            global lp_b
            
            dim = 3; % 3D
            uff_t = obj.u0.eval(t);      % feedforward control
            nq = obj.nq;
            nv = obj.nv;            
            nc = 4; %obj.robot.getNumContactPairs();
            qa = obj.qa; 
            
            % desired actuator state
            x_des_t = obj.x_des.eval(t);
            q_des_t = x_des_t(1:nq);
            qd_des_t = x_des_t(nq+(1:nv));       
            
            % current actuator state
            q_t = x(1:nq);
            qd_t = x(nq+(1:nv));
            
            switch obj.CTYPE
                
                case 'legpd'
                    
                    % desired leg state in body frame
                    qleg_des_t = zeros(dim,  nc);
                    vleg_des_t = qleg_des_t;
                    
                    fkopt.base_or_frame_id = obj.robot.findLinkId('Chassis');
                    kinsol_des = obj.robot.doKinematics(q_des_t, qd_des_t);
                    for i = 1:nc
                        [qleg_des_t(:, i), Jleg_des] = obj.robot.forwardKin(kinsol_des, ...
                            obj.robot.findLinkId(obj.robot.legs{i}), lp_b(i,:)', fkopt);
                        vleg_des_t(:,i) = Jleg_des*qd_des_t;
                    end
                    
                    % current leg state in body frame
                    qleg_t = zeros(dim, nc);
                    vleg_t = qleg_t;
                    
                    kinsol = obj.robot.doKinematics(q_t, qd_t);
                    for i = 1:nc
                        [qleg_t(:, i), JlegB] = obj.robot.forwardKin(kinsol, ...
                            obj.robot.findLinkId(obj.robot.legs{i}), lp_b(i,:)', fkopt);
                        vleg_t(:,i) = JlegB*qd_t;
                        
                    end
                    
                    % leg state error
                    q_err_t = qleg_des_t([1,3], :) - qleg_t([1,3], :);
                    v_err_t = vleg_des_t([1,3], :) - vleg_t([1,3], :);
                    
                    % gain matrix
                    K_t = obj.K;
                    K_t([1, 6, 7, 8], :) = -K_t([1, 6, 7, 8], :);
                    
                    % control 
                    u = uff_t + K_t*[q_err_t(:);v_err_t(:)];
                    
%                     fprintf('Leg Error: %f at %f sec \r', max(abs(q_err_t(:))), t)                   
                    
                case 'actlqr'
                    
                    q_err_t = q_des_t(qa) - q_t(qa);
                    qd_err_t = qd_des_t(qa) - qd_t(qa);  
                    x_err_t = [q_err_t, qd_err_t]';
                    K_t = obj.K;                     
                    
                    u = uff_t + K_t*x_err_t(:);
                    
%                     q_err_t(1)^2*100
%                     qd_err_t(1)^2*100
                    
%                     fprintf('\n Act Error: %f mm, %f mm/ms at %f sec \r', max(abs(q_err_t(:))), max(abs(qd_err_t(:))), t)
%                     fprintf('\n Voltage: %f V, %f mm/ms at %f sec \r',  max(u)^2, t)
                    
            end
            
            
            % Input Limits
            u(u < obj.act.dummy_bender(1).dp.Vg) = obj.act.dummy_bender(1).dp.Vg;
            u(u > obj.act.dummy_bender(1).dp.Vb) = obj.act.dummy_bender(1).dp.Vb; 
            
        end
        
        
    end
    
end