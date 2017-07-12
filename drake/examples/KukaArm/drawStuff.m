load SampleTrajectories

options.with_weight = true;
options.with_shelf_and_boxes = true;
p = KukaArm(options);

N = length(tsamp1);
qsamp1 = xsamp1(1:7,:);
qtraj1 = PPTrajectory(foh(tsamp1,qsamp1));
qsamp2 = xsamp2(1:7,:);
qtraj2 = PPTrajectory(foh(tsamp2,qsamp2));

for k = 1:N
    ks1 = p.doKinematics(qsamp1(:,k));
    hand1(:,k) = p.forwardKin(ks1,findLinkId(p,p.hand_name),[0;0;0]);
    ks2 = p.doKinematics(qsamp2(:,k));
    hand2(:,k) = p.forwardKin(ks2,findLinkId(p,p.hand_name),[0;0;0]);
end

v = p.constructVisualizer();
qtraj1 = qtraj1.setOutputFrame(v.getInputFrame());
qtraj2 = qtraj2.setOutputFrame(v.getInputFrame());
gl = drake.util.BotLCMGLClient(lcm.lcm.LCM.getSingleton(),'kuka');

gl.glLineWidth(6.5);
gl.glColor3f(0,0,1);
gl.plot3(hand1(1,:),hand1(2,:),hand1(3,:));
gl.glLineWidth(6.5);
gl.glColor3f(1,0,0);
gl.plot3(hand2(1,:),hand2(2,:),hand2(3,:));
gl.switchBuffers();

v.playback(qtraj1,struct('slider',true));
v.playback(qtraj2,struct('slider',true));