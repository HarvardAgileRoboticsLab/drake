classdef FoamyPendulumVisualizer < Visualizer
    % Implements the draw function for the McFoamy RC airplane
    
    methods
        function obj = FoamyPendulumVisualizer(plant)
          typecheck(plant,'FoamyPendulumPlant');
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
            
          obj.model = plant.plane_rbm;
          %obj.geom = RigidBodyMesh('flat-plate.obj', [0 0.44 0.2]', [0 0 0]');
          %obj.geom.scale = 0.001;
          %obj.geom.c = [0.6,0,0];
          %obj.geom(2) = RigidBodyMesh('sphere.obj', [0 0 0]', [0 0 0]');
          %obj.geom(2).scale = 1;
          %obj.geom(2).c = [0,0,0];
          %obj.geom = {RigidBodyMesh('flat-plate.obj', [0 0.44 0.2]', [0 0 0]')...
          %             RigidBodyMesh('sphere.obj', [0 0 0]', [0 0 0]') };
          obj.geom = {obj.model.body(2).visual_geometry,...
              obj.model.body(3).visual_geometry,...
              obj.model.body(4).visual_geometry};
          %obj.geom{1}.scale = 0.001;
          %obj.geom{1}.c = [0.6,0,0];
          %obj.geom{2}.c = [0,0,0];
          obj = update(obj,plant);
          manip = plant.plane_rbm;
         obj.draw_msg = drake.lcmt_viewer_draw();
         nb = getNumBodies(manip);
         obj.draw_msg.num_links = nb;
         obj.draw_msg.link_name = {manip.body.linkname};
         obj.draw_msg.robot_num = [manip.body.robotnum];
         obj.draw_msg.position = single(zeros(nb,3));
         obj.draw_msg.quaternion = single(zeros(nb,4));

          draw(obj,0,getZeroConfiguration(plant));
        end
  
  function obj = update(obj,plant)
      use_collision_geometry = 0;
      manip = plant.plane_rbm;
%      obj = update@Visualizer(obj,plant.getOutputFrame);
      
      lc = lcm.lcm.LCM.getSingleton();
      vr = drake.lcmt_viewer_load_robot();
      vr.num_links = getNumBodies(manip);
      vr.link = javaArray('drake.lcmt_viewer_link_data',vr.num_links);
      for i=1:vr.num_links
        b = getBody(manip,i);
        link = drake.lcmt_viewer_link_data();
        link.name = b.linkname;
        link.robot_num = b.robotnum;
        if use_collision_geometry
          link.num_geom = length(b.collision_geometry);
        else
          link.num_geom = length(b.visual_geometry);
        end
        if (link.num_geom>0)
          link.geom = javaArray('drake.lcmt_viewer_geometry_data',link.num_geom);
        end
        for j=1:link.num_geom
          if use_collision_geometry
            link.geom(j) = serializeToLCM(b.collision_geometry{j});
          else
            link.geom(j) = serializeToLCM(b.visual_geometry{j});
          end
        end
        vr.link(i) = link;
      end
      
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
      
      %obj.draw_msg.position = q(1:3)';
      %obj.draw_msg.quaternion = q(4:7)';
      
        obj.draw_msg.position(2,:) =  q(1:3)';
        obj.draw_msg.quaternion(2,:) = q(4:7)';
        obj.draw_msg.position(3,:) = q(8:10)';
        obj.draw_msg.quaternion(3,:) = q(11:14)';
        obj.draw_msg.position(4,:) = q(8:10)';
        obj.draw_msg.quaternion(4,:) = q(11:14)';
      
      lc = lcm.lcm.LCM.getSingleton();
      lc.publish('DRAKE_VIEWER_DRAW',obj.draw_msg);
      
    end
    end
    
  properties
    geom;
    viewer_id;
    draw_msg;
    status_agg;
    model;
  end
    
end
