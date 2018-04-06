function plot_traj(veh,ts,vec_length)
%Veh= vehicla struct output from traj_from_log
%ts = [t0 t1] Indexes between which to plot
%vec_length = length of vector for drawing coordinate axes for (preference)

%This function will take simulation/flight log data and plot between
%desired points the measured/desired positions and orientations.


figure(11)

plot_des = 1;
plot_meas = 1;
plot_origin = 0;
plot_heading = 0;

%x_init = veh.x_init(ts(1):ts(2),:);
%xm = veh.x_meas(ts(1):ts(2),:);

%xm = [1*veh.x(ts(1):ts(2),:)';1*veh.y(ts(1):ts(2),:)';1*veh.z(ts(1):ts(2),:)';...
%    -1*veh.q1(ts(1):ts(2),:)';1*veh.q0(ts(1):ts(2),:)';-1*veh.q3(ts(1):ts(2),:)';1*veh.q2(ts(1):ts(2),:)']';
%xd(1:end,3) = xd(1:end,3) + x_init(1:end,3);
%xd(1:end,1:2) = xd(1:end,1:2) - x_init(1:end,1:2);
%xd(1:end,3) = xd(1:end,3)*-1;

if(plot_origin)
    
    plot3([0,.05],[0 0],[0 0],'r-','linewidth',2);
    hold on
    plot3([0,0],[0 -.05],[0 0],'g-','linewidth',2);
    plot3([0,0],[0 0],[0 -.05],'b-','linewidth',2);
end

if(plot_des)
    xd = veh.xd(ts(1):ts(2),:);
    %xd(1:end,3) = xd(1:end,3)*-1;
    plot3(xd(:,1),xd(:,2),xd(:,3),'k-')
    hold on
    plot3(xd(1,1),xd(1,2),xd(1,3),'go')
    plot3(xd(end,1),xd(end,2),xd(end,3),'ro')
    axis equal
end

if(plot_meas)
    xm = veh.x_meas(ts(1):ts(2),:);
    plot3(xm(:,1),xm(:,2),xm(:,3),'c-')
    hold on
    plot3(xm(1,1),xm(1,2),xm(1,3),'go')
    plot3(xm(end,1),xm(end,2),xm(end,3),'ro')
    axis equal
end

[N,~] = size(xm);

for i = 1:N

    
    if(plot_des)

        q_res = xd(i,4:7)';

        x_dir = qrotate(q_res,[vec_length;0;0]);
        y_dir = qrotate(q_res,[0;vec_length;0]);
        z_dir = qrotate(q_res,[0;0;vec_length]);

        plot3([xd(i,1),xd(i,1)+x_dir(1)],[xd(i,2),xd(i,2)+x_dir(2)],[xd(i,3),xd(i,3)+x_dir(3)],'r-');
        plot3([xd(i,1),xd(i,1)+y_dir(1)],[xd(i,2),xd(i,2)+y_dir(2)],[xd(i,3),xd(i,3)+y_dir(3)],'g-');
        plot3([xd(i,1),xd(i,1)+z_dir(1)],[xd(i,2),xd(i,2)+z_dir(2)],[xd(i,3),xd(i,3)+z_dir(3)],'b-');



        if(plot_heading)
             plot3([xm(i,1),xm(i,1)+vec_length*cos(x_init(i,4))],[xm(i,2),xm(i,2)+vec_length*sin(x_init(i,4))],[xd(i,3),xd(i,3)],'k.-');
        end

    end
    if(plot_meas) 

        q_resm = xm(i,4:7)';

        xm_dir = qrotate(q_resm,[vec_length;0;0]); 
        ym_dir = qrotate(q_resm,[0;vec_length;0]);
        zm_dir = qrotate(q_resm,[0;0;vec_length]);

        plot3([xm(i,1),xm(i,1)+xm_dir(1)],[xm(i,2),xm(i,2)+xm_dir(2)],[xm(i,3),xm(i,3)+xm_dir(3)],'r-');
        plot3([xm(i,1),xm(i,1)+ym_dir(1)],[xm(i,2),xm(i,2)+ym_dir(2)],[xm(i,3),xm(i,3)+ym_dir(3)],'g-');
        plot3([xm(i,1),xm(i,1)+zm_dir(1)],[xm(i,2),xm(i,2)+zm_dir(2)],[xm(i,3),xm(i,3)+zm_dir(3)],'b-');


    end

end



hold off
end

function rrot = qrotate(q,r)
    %Rotate vector r by quaternion q
    rrot = r + 2*cross(q(2:4),(cross(q(2:4),r) + q(1)*r));
end

function q_res = qrotate_eul(q,r,theta)

qrot = [cos(theta/2);sin(theta/2)*r(1);sin(theta/2)*r(2);sin(theta/2)*r(3)];
qrot = [qrot(1);qrot(2:4)];
q_res = qmultiply(qrot,q);


end

function q = qmultiply(q1,q2)

    %Quaternion multiplication
    s1 = q1(1);
    v1 = q1(2:4);
    s2 = q2(1);
    v2 = q2(2:4);
    q = [s1*s2 - v1'*v2; s1*v2 + s2*v1 + cross(v1,v2)];

end