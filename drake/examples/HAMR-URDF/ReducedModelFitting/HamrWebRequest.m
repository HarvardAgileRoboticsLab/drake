classdef HamrWebRequest < DrakeSystem
    
    
    properties (SetAccess=private)
        
        url;
        media_type;
        
        nq;
        nv;
        nu;
        ihip;
        perm; 
        
    end
    
    
    methods
        function obj = HamrWebRequest(rwact,r)
            
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
            
        end
        
        
        function y = output(obj, t, ~, x)
            disp(t)
            
            xhip = x(obj.ihip);
            
            temp = struct( ...
                'q_pos', x(1:6), ...
                'q_hip', xhip(obj.perm), ...
                'qdot_pos', x(obj.nq + (1:6)), ...
                'qdot_hip', xhip(obj.nu + obj.perm) ...
                );
            
            %             disp('Web request...')
            y_josh = webwrite(obj.url, temp, obj.media_type);
            y = y_josh([1, 5, 2, 6, 3, 7, 4, 8]);         
            
        end
        
        
    end
    
end