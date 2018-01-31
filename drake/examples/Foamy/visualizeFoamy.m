function [] =  visualizeFoamy(obj,xtraj,slider)

type = class(obj);

if (type == "FoamyPlant")
       v = FoamyVisualizer(obj);
       v.playback(xtraj,struct('slider',slider));
    
elseif (type == "FoamyPendulumPlant")
       v = FoamyPendulumVisualizer(obj);
       v.playback(xtraj,struct('slider',slider));
    
elseif (type == "HybridFoamyPlant")
            xtraj_mode_1 = xtraj.traj{1}.trajs{2};
            xtraj_mode_1=xtraj_mode_1.setOutputFrame(obj.getOutputFrame);
            xtraj_mode_2 = xtraj.traj{2}.trajs{2};
            xtraj_mode_2=xtraj_mode_2.setOutputFrame(obj.getOutputFrame);
            xtraj = xtraj_mode_1.append(xtraj_mode_2);
            v = FoamyVisualizer(obj);
            v.playback(xtraj,struct('slider',slider));
    
elseif (type == "HybridFoamyPendulumPlant")
            xtraj_mode_1 = xtraj.traj{1}.trajs{2};
            xtraj_mode_1=xtraj_mode_1.setOutputFrame(obj.getOutputFrame);
            xtraj_mode_2 = xtraj.traj{2}.trajs{2};
            xtraj_mode_2=xtraj_mode_2.setOutputFrame(obj.getOutputFrame);
            xtraj = xtraj_mode_1.append(xtraj_mode_2);
            v = FoamyPendulumVisualizer(obj);
            v.playback(xtraj,struct('slider',slider));
    
else
    disp("Please call the function for the correct class");
    

        


end

