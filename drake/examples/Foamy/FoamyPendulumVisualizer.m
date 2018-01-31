classdef FoamyPendulumVisualizer < Visualizer
    % Implements the draw function for the McFoamy RC airplane
    
    methods
        function obj = FoamyPendulumVisualizer(plant)
            
            %typecheck(plant,'FoamyPendulumPlant');
            checkDependency('lcm');
            
            obj = obj@Visualizer(plant.getOutputFrame);
            
            lc = lcm.lcm.LCM.getSingleton();
            obj.status_agg = lcm.lcm.MessageAggregator();
            lc.subscribe('DRAKE_VIEWER_STATUS',obj.status_agg);
            
            % check if there is a viewer already running
            [~,ck] = system('ps ax 2> /dev/null | grep -i "drake-visualizer" | grep -c -v grep');
            if (str2num(ck)<1)
                % try launching director first
                if exist(fullfile(drake_get_bin_path,'drake-visualizer'))
                    disp('attempting to launch the drake director')
                    retval = systemWCMakeEnv([fullfile(drake_get_bin_path,'drake-visualizer'),' &> drake-visualizer.out &']);
                    
                    if isempty(obj.status_agg.getNextMessage(10000)) % wait for viewer to come up
                        type drake-visualizer.out
                        error('FoamyVisualizer:AutostartFailed','Failed to automatically start up a viewer (or to receive the ack, see https://github.com/RobotLocomotion/drake/issues/317)');
                    end
                end
            end
            
            obj.model = plant;
            obj.plane_geom = RigidBodyMesh('flat-plate.obj', [-0.7 0.4375 0.1]', [0 0 0]');
            obj.plane_geom.scale = 0.001;
            obj.plane_geom.c = [0.6,0,0];

            obj.pendulum_geom = RigidBodyCylinder(.01, .6);
            obj.pendulum_geom.c = [0,0,0.8];
            
            obj = update(obj,plant);
            
            obj.draw_msg = drake.lcmt_viewer_draw();
            obj.draw_msg.num_links = 2;
            obj.draw_msg.link_name = {'foamy','pendulum'};
            obj.draw_msg.robot_num = [1 1];
            obj.draw_msg.position = single(zeros(2,3));
            obj.draw_msg.quaternion = single(zeros(2,4));
            
            draw(obj,0,getZeroConfiguration(plant));
        end
        
        function obj = update(obj,plant)
            
            lc = lcm.lcm.LCM.getSingleton();
            vr = drake.lcmt_viewer_load_robot();
            vr.num_links = 2;
            vr.link = javaArray('drake.lcmt_viewer_link_data',vr.num_links);
            
            link1 = drake.lcmt_viewer_link_data();
            link1.name = 'foamy';
            link1.robot_num = 1;
            link1.num_geom = 1;
            link1.geom = javaArray('drake.lcmt_viewer_geometry_data',1);
            link1.geom(1) = serializeToLCM(obj.plane_geom);
            vr.link(1) = link1;
            
            link2 = drake.lcmt_viewer_link_data();
            link2.name = 'pendulum';
            link2.robot_num = 1;
            link2.num_geom = 1;
            link2.geom = javaArray('drake.lcmt_viewer_geometry_data',1);
            link2.geom(1) = serializeToLCM(obj.pendulum_geom);
            vr.link(2) = link2;
            
            lc.publish('DRAKE_VIEWER_LOAD_ROBOT',vr);
            
            if (false) % the message aggregator is missing valid acks
                % listen for acknowledgement
                ack = obj.status_agg.getNextMessage(5000);
                obj.status_agg.numMessagesAvailable()
                if isempty(ack)
                    error('FoamyVisualizer:LoadRobotFailed','Did not receive ack from viewer');
                else
                    msg = drake.lcmt_viewer_command(ack.data);
                    
                end
            end
        end
        
        function drawWrapper(obj,t,y)
            draw(obj,t,y);
        end
        
        function draw(obj,t,q)
            obj.draw_msg.timestamp = int64(t*1000000);
            
            obj.draw_msg.position = [q(1:3)'; q(8:10)'];
            obj.draw_msg.quaternion = [q(4:7)'; q(11:14)'];
            
            lc = lcm.lcm.LCM.getSingleton();
            lc.publish('DRAKE_VIEWER_DRAW',obj.draw_msg);
            
        end
    end
    
    properties
        plane_geom;
        pendulum_geom;
        viewer_id;
        draw_msg;
        status_agg;
        model;
    end
    
end
