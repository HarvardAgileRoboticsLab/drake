%Initial State, Controls, and Sensors
x0 = [0 0 0 1/sqrt(2) 0 0 1/sqrt(2) 0 0 0 0 0 0]';
u0 = [0 127 127 127]';
y0 = foamy_sensors(x0,u0);

%Listen for HIL_ACTUATOR_CONTROLS messages
% receiver = dsp.UDPReceiver('LocalIPPort',14560);
% setup(receiver);

%Send HIL_SENSOR, HIL_STATE_QUATERNION, and HIL_GPS messages
sender = dsp.UDPSender('RemoteIPPort',14560);

while true
    ydata = mavlink_pack_sensors_mex(y0);
    step(sender, ydata);
    xdata = mavlink_pack_state_mex(x0,y0);
    step(sender, xdata);
    pause(.002);
%     packet = step(receiver);
%     if ~isempty(packet)
%         [u,t] = mavlink_parse_control_mex(packet);
%     end
end

