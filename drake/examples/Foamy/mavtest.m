%--- Simulation Parameters ---%
ts = .02; %Simulation timestep in seconds
x = [0 0 0 0 1/sqrt(2) 1/sqrt(2) 0 0 0 0 0 0 0]'; %Initial State (pointing North)
u = [0 0 0 0]'; %Initial controls

plant = FoamyPlant();
vis = FoamyVisualizer(plant);
vis.draw(0,x);

%Listen for HIL_ACTUATOR_CONTROLS messages
receiver = dsp.UDPReceiver('LocalIPPort',54321);
setup(receiver);

%Send HIL_SENSOR, HIL_STATE_QUATERNION, and HIL_GPS messages
sender = dsp.UDPSender('RemoteIPPort',14560,'LocalIPPortSource','Property','LocalIPPort',54321);

while true
    tic
    
    %Simulate forward one timestep with midpoint integration
    [x, xdot] = foamy_rk4(x,u,ts);
    
    %Calculate sensor measurements and add noise
    %(PX4 doesn't work without some noise)
    y = foamy_sensors(x,xdot)+0.0001*randn(17,1);

    %Communicate over MAVLink
    ydata = mavlink_pack_sensors_mex(y);
    step(sender, ydata);
    xdata = mavlink_pack_state_mex(x,y);
    step(sender, xdata);
    gdata = mavlink_pack_gps_mex(y);
    step(sender, gdata);
    packet = step(receiver);
    if ~isempty(packet)
        if packet(8) == 93
            u = mavlink_parse_control_mex(packet);
        end
    end
    
    vis.draw(0,x);
    
    %Wait so we stay synchronized
    dt = toc;
    pause(ts-dt);
end

