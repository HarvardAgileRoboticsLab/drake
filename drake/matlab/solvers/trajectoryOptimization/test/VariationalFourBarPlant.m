classdef VariationalFourBarPlant< RigidBodyManipulator
    
    properties (SetAccess = protected, GetAccess = public)
        x0
        nKL = 1;  % number of kinematic loops
    end
    
    methods
        
        function obj=VariationalFourBarPlant(urdf, options)
            
            typecheck(urdf,'char');
            %             options.floating = false;
            %             options.use_bullet = false;
            
            obj = obj@RigidBodyManipulator(urdf,options);
            nQ = obj.getNumPositions();
            q0 = ones(nQ, 1)*pi/2;
            obj.x0 = [q0; 0*q0];
            
            
            klvars = nQ;
            
            cnstr_opts.grad_level = 2;
            cnstr_opts.grad_method = 'user';
            
            klconst = FunctionHandleConstraint(6*zeros(obj.nKL, 1), ...
                6*zeros(obj.nKL, 1), klvars, @obj.kl_const_fun, cnstr_opts);
            obj = obj.addPositionEqualityConstraint(klconst, 1:4);  
            
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
            
            nQ = obj.getNumPositions();
            kinsol = doKinematics(obj,q,zeros(nQ,1), struct('compute_gradients', true));            
            
            body1 = 1;
            body2 = 5;
            [x1, J1, dJ1] = obj.forwardKin(kinsol, body1, zeros(3,1), struct('rotation_type', 1));
            [x2, J2, dJ2] = obj.forwardKin(kinsol, body2, zeros(3,1), struct('rotation_type', 1));
%             [xrel, J, dJ] = obj.forwardKin(kinsol, body1, zeros(3,1), kinopt); 
            
%             s = xrel(4); v = xrel(5:7); 
%             r = v'*v/s^2; 
%             tht = 2*atan2(sqrt(v'*v), sqrt(s^2)); 
%             axis = v/sin(tht/2); 
%             rad2deg(tht)
%             f = x(4) - tht; 
%             df = [0, 0, 0, 1] - (2/s)*(1/(sqrt(r)*(1+r)))*[-r, v'/s]*J(4:7, :); 
%             df2 = -(2/s)*(1/(sqrt(r)*(1+r)))*[-r, v'/s]*dJ(4:7,:);
%             f = [q(4) - q(2); q(3) - q(1)];  
%             df = [0, -1, 0, 1; -1, 0, 1, 0]; 
%             df2 = zeros(2, nQ*nQ); 
%             

        %             
            f = x2 - x1;
            df = J2 - J1;
            df2 = dJ2 - dJ1;
            
            
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
            x0 = obj.x0;
        end
    end
    
    
    
end
