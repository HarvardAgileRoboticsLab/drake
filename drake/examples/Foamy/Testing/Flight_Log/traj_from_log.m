function [xtraj,utraj] = traj_from_log(log_name,p)

%log name is a character array of the log name
%p is a foamy plant object


csv_files = dir('*.csv');

dir_log = strcat('flight_logs_csv/',log_name);
dir_log2 = strcat('./',dir_log);

if(isempty(csv_files) && ~exist(dir_log2,'dir'))
    ulog_file_name = log_name + ".ulg";
    system('ulog2csv '+ ulog_file_name);
    mkdir(dir_log);
    temp = dir('*.csv');
    N = length(temp);
    
   for i = 1:N
       movefile(temp(i).name, dir_log)
   end
   % mkdir('flight_logs_csv/'+log_name);
end

%cd (dir_log)

local_position_filename = log_name + "_vehicle_local_position_0.csv";
vehicle_attitude_filename = log_name + "_vehicle_attitude_0.csv";
actuator_controls_filename = log_name + "_actuator_controls_0_0.csv";

local_pos = csvread(strcat(dir_log2,"/"+local_position_filename),1,0);
veh_att = csvread(strcat(dir_log2,"/"+vehicle_attitude_filename),1,0);
controls = csvread(strcat(dir_log2,"/"+actuator_controls_filename),1,0);

size(local_pos(:,6))

vehicle = {};
vehicle.x = local_pos(:,6);
vehicle.y = local_pos(:,7);
vehicle.z = -1*local_pos(:,8);

vehicle.vx = local_pos(:,12);
vehicle.vy = local_pos(:,13);
vehicle.vz = local_pos(:,14);

vehicle.q0 = veh_att(:,5);
vehicle.q1 = veh_att(:,6);
vehicle.q2 = veh_att(:,7);
vehicle.q3 = veh_att(:,8);

vehicle.w_roll = veh_att(:,2);
vehicle.w_pitch = veh_att(:,3);
vehicle.w_yaw = veh_att(:,4);

initial_time = veh_att(1,1);
att_time = 10^(-6)*(veh_att(:,1)-veh_att(1,1));
pos_time = 10^(-6)*(local_pos(:,1)-veh_att(1,1));

vehicle.x = (vehicle.x-vehicle.x(1));
vehicle.y = -1*(vehicle.y-vehicle.y(1));

%xtraj = [local_pos(:,1),vehicle.x,vehicle.y,vehicle.z,...
%    vehicle.q0,vehicle.q1,vehicle.q2,vehicle.q3,...
%    vehicle.vx,vehicle.vy,vehicle.vz,vehicle.w_roll,vehicle.w_pitch,vehicle.w_yaw];


xtraj_val = [vehicle.x,vehicle.y,vehicle.z,...
    vehicle.vx,vehicle.vy,vehicle.vz];
xtraj_ts = [local_pos(:,1)];

x_pp = foh(pos_time',xtraj_val');

x_vals = ppval(x_pp,att_time);

%qs = [vehicle.q0,vehicle.q1,vehicle.q2,vehicle.q3];

%qs = qmultiply(qs,[0,1,0,0]);

%xtraj = [x_vals(1:3,:)',vehicle.q0,vehicle.q1,vehicle.q2,vehicle.q3,...
%    x_vals(4:6,:)',vehicle.w_roll,vehicle.w_pitch,vehicle.w_yaw];

xtraj = [x_vals(1:3,:)',-1*vehicle.q1,vehicle.q0,-1*vehicle.q3,-1*vehicle.q2,...
    x_vals(4:6,:)',vehicle.w_roll,vehicle.w_pitch,vehicle.w_yaw];

ts = att_time;

pp = foh(ts',xtraj');

xtraj = PPTrajectory(pp);
%xtraj = xtraj.setOutputFrame(p.getOutputFrame);
xtraj = xtraj.setOutputFrame(p.getStateFrame);



%figure(1)
%plot3(vehicle.x,vehicle.y,vehicle.z)


%xtraj = 0;
utraj = 0;

end

function q = qmultiply(q1, q2)
    %Quaternion multiplication
    s1 = q1(1);
    v1 = q1(2:4);
    s2 = q2(1);
    v2 = q2(2:4);
    q = [s1*s2 - v1'*v2; s1*v2 + s2*v1 + cross(v1,v2)];
end

function qc = qconj(q)
    %Quaternion conjugate
    qc = [q(1); -q(2:4)];
end

