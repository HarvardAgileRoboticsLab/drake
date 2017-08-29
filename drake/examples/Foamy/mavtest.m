%--- Simulation Parameters ---%
ts = .02; %Simulation timestep in seconds
x = [0 0 0 0 1/sqrt(2) 1/sqrt(2) 0 0 0 0 0 0 0]'; %Initial State (pointing North)
u = [0 0 0 0]'; %Initial controls

%Send HIL_SENSOR, HIL_STATE_QUATERNION, and HIL_GPS messages
sender = dsp.UDPSender('RemoteIPPort',14560,'LocalIPPortSource','Property','LocalIPPort',54321);

%Listen for HIL_ACTUATOR_CONTROLS messages
receiver = dsp.UDPReceiver('LocalIPPort',54321);
setup(receiver);

while true
    tic
    
    %Simulate forward one timestep with Euler integration
    xdot = foamy_dynamics_mex(0,x,u);
    x = x + ts*xdot;
    
    %Make sure we don't go through the ground
    x(3) = max(0,x(3));
    if(x(3) == 0)
        x(10) = max(0,x(10));
        xdot(3) = max(0,xdot(3));
        xdot(10) = max(0,xdot(10));
    end
    
    %Make sure quaternion is normalized
    x(4:7) = x(4:7)/sqrt(x(4:7)'*x(4:7));
    
    %Calculate sensor measurements and add noise
    y = foamy_sensors(x,xdot)+0.001*randn(17,1);

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
    
    %Wait so we stay synchronized
    dt = toc;
    pause(ts-dt);
end

