%--- Simulation Parameters ---%
ts = .02; %Simulation timestep in seconds
x = [0 0 0 0 1 0 0 0 0 0 0 0 0]'; %Initial State (pointing East)
u = [0 0 0 0]'; %Initial controls

plant = FoamyPlant();
vis = FoamyVisualizer(plant);
vis.draw(0,x);

%Listen for HIL_ACTUATOR_CONTROLS messages
receiver = dsp.UDPReceiver('LocalIPPort',54321);

%Send HIL_SENSOR, HIL_STATE_QUATERNION, and HIL_GPS messages
sender = dsp.UDPSender('RemoteIPPort',14560,'LocalIPPortSource','Property','LocalIPPort',54321);

sim_timer = tic;
gps_timer = tic;
loop_timer = tic;
while true
    
    %Try to receive control inputs over MAVLink
    packet = step(receiver);
    if ~isempty(packet)
        if packet(8) == 93
            u = mavlink_parse_control_mex(packet);
        end
    end
    
    %Simulate forward one timestep
    dt = toc(loop_timer);
    [x, xdot] = foamy_midpoint(x,u,dt);
    loop_timer = tic;
    
    %Calculate sensor measurements and add noise
    %(PX4 doesn't work without some noise)
    y = foamy_sensors(x,xdot);
    yn = y + [1e-7; 1e-7; .01; 1e-5*ones(14,1)].*randn(17,1);
    
    %Send updated state + sensor measurements over MAVLink
    now = toc(sim_timer);
    ydata = mavlink_pack_sensors_mex(yn,now);
    step(sender, ydata);
    xdata = mavlink_pack_state_mex(x,y,now);
    step(sender, xdata);
    if (toc(gps_timer) > 0.1)
        gdata = mavlink_pack_gps_mex(yn,now);
        step(sender, gdata);
        gps_timer = tic;
    end
%     vdata = mavlink_pack_vision_mex(x);
%     step(sender, vdata);
    
    vis.draw(0,x);
end

