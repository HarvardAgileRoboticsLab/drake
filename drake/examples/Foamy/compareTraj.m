function compareTraj(xtraj1,xtraj2,utraj1,utraj2)
tspan = min(xtraj1.tspan(2),xtraj2.tspan(2));
ts = linspace(0,tspan);


x1 = xtraj1.eval(ts);
quats1 = x1(5:8,:);
eul1 = quat2eul(quats1');
xs_1 = [x1(2:4,:);eul1';x1(9:end,:)];
us_1 = utraj1.eval(ts);

x2 = xtraj2.eval(ts);
quats2 = x2(4:7,:);
eul2 = quat2eul(quats2');
xs_2 = [x2(1:3,:);eul2';x2(8:end,:)];
us_2 = utraj2.eval(ts);

[m_x,~] = size(xs_2);

[m_u,~] = size(us_2);

figure(2);

titles = {'x','y','z','th1','th2','th3','dx','dy','dz','w1','w2','w3',...
    'u1','u2','u3','u4'};

for i = 1:m_x

    subplot(3,6,i)
    %plot(ts,xs_1(i,:));
    hold on
    plot(ts,xs_2(i,:));
    title(titles(i))
    hold off
end

for i = 1:m_u
    subplot(3,6,i+m_x)
    %plot(ts,us_1(i,:));
    hold on
    plot(ts,us_2(i,:));
    title(titles(m_x+i))
    hold off
end

end

