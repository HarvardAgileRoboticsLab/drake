%Listen for HIL_ACTUATOR_CONTROLS messages
receiver = dsp.UDPReceiver('LocalIPPort',14560);
setup(receiver);

%Send HIL_SENSOR, HIL_STATE_QUATERNION, and HIL_GPS messages
sender = dsp.UDPSender('RemoteIPPort',14560);

while true
    packet = step(receiver);
    if ~isempty(packet)
        [u,t] = mavlink_parse_control_mex(packet);
    end
end

