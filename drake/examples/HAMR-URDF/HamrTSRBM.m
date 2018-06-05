classdef HamrTSRBM < TimeSteppingRigidBodyManipulator
    
    properties (SetAccess = protected, GetAccess = public)
        HamrRBM = [];
        %         q0
        %         grav = [0; 0; -9.81e-3];
        %         valid_loops;
        %         legs = {'FLL4', 'RLL4', 'FRL4', 'RRL4'};
        %         ULIM = 0.3;
        %         stiffness_multipliers = ones(4,1);      % default (lift act, swing act, lift flex, swing flex);
    end
    
    
    methods
        
        function obj=HamrTSRBM(manip, dt, options)
            
            typecheck(manip,'RigidBodyManipulator');
            obj = obj@TimeSteppingRigidBodyManipulator(manip,dt,options);
            obj.HamrRBM = manip;      
            
        end
        
        function [y, Jy] = output(obj, t, x, u)     % call RBM's output function
            [y, Jy] = obj.HamrRBM.output(t, x, u);
        end
                
        function nActuatedDOF = getNumActuatedDOF(obj)
            nActuatedDOF = numel(obj.getActuatedJoints());
        end
        
        function obj = setInitialState(obj,x0)
            typecheck(x0,'double');
            sizecheck(x0,obj.getNumStates());
            obj.HamrRBM = obj.HamrRBM.setInitialState(x0);
        end
        
        function x0 = getInitialState(obj)
            x0 = obj.HamrRBM.x0;
        end
        
    end
end

