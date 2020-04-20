classdef HamrWebRequest < DrakeSystem
    
    properties (SetAccess=private)
        
        r;
        
        url;
        media_type;
        opt_2x;
        
        nq;
        nv;
        nu;
        ihip;
        perm; 
        
    end
    
    
    methods
        function obj = HamrWebRequest(rwact,r, opt_2x)
            
            obj = obj@DrakeSystem(0,0,rwact.getNumOutputs, ...
                rwact.getNumInputs,true,true);
            
            obj = setInputFrame(obj,rwact.getOutputFrame());
            obj = setOutputFrame(obj,rwact.getInputFrame());

            obj.nq = r.getNumPositions();
            obj.nv = r.getNumVelocities();
            obj.nu = r.getNumInputs();
            obj.ihip = find(contains(r.getStateFrame().getCoordinateNames(), ...
                {'l2', 's2'}) == 1);  % index of SFB inputs
            obj.perm = [1, 3, 5, 7, 2, 4, 6, 8];
            
            obj.url = 'http://35.245.168.29:8000';
            obj.media_type = weboptions('MediaType', 'applications/json');
            
            obj.opt_2x = opt_2x;
            
            obj.r = r; 
            
        end
        
        
        function y = output(obj, t, ~, x)
            persistent ytm1

            
            if obj.opt_2x 
                flag = mod(t, 2*obj.r.timestep); 
            else
                flag = 0;
            end
            
            xhip = x(obj.ihip);
            
            temp = struct( ...
                'q_pos', [-x(1); x(2:6)], ...
                'q_hip', xhip(obj.perm), ...
                'qdot_pos', [-x(obj.nq + 1); x(obj.nq + (2:6))], ...
                'qdot_hip', xhip(obj.nu + obj.perm) ...
                );

            if flag == 0  
                disp(t)                       
                y_josh = webwrite(obj.url, temp, obj.media_type);
                y = y_josh([1, 5, 2, 6, 3, 7, 4, 8]);              
                ytm1 = y; 
%                 obj = obj; 
%                 disp(y)
            
            else
                y = ytm1;
%                 disp(y)
            end
            
            
        end
        
        
    end
    
end