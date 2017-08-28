%Initial State, Controls, and Sensors
x0 = [0 0 0 1/sqrt(2) 0 0 1/sqrt(2) 0 0 0 0 0 0]';
u0 = [0 127 127 127]';
y0 = foamy_sensors(x0,u0);

%Send HIL_SENSOR, HIL_STATE_QUATERNION, and HIL_GPS messages
sender = dsp.UDPSender('RemoteIPPort',14560);

%Listen for HIL_ACTUATOR_CONTROLS messages
% receiver = dsp.UDPReceiver('LocalIPPort',14560);
% setup(receiver);

while true
    y = y0+0.001*randn(15,1);
    x = x0+0.001*randn(13,1);
    ydata = mavlink_pack_sensors_mex(y);
    step(sender, ydata);
    xdata = mavlink_pack_state_mex(x,y);
    step(sender, xdata);
    gdata = mavlink_pack_gps_mex(y);
    step(sender, gdata);
    vdata = mavlink_pack_vision_mex(x);
    step(sender, vdata);
    pause(.004);
%     packet = step(receiver);
%     if ~isempty(packet)
%         [u,t] = mavlink_parse_control_mex(packet);
%     end
end

