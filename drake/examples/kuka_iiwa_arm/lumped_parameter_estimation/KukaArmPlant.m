classdef KukaArmPlant < Manipulator 
  
  properties
    % parameters from KUKA CAD Model 
    
    %(except inertias are now relative to the
    % joints)
    % axis)%[To Be Checked]
    
    l1x = 0, l1y = 0, l1z = 0.1575;
    l2x = 0, l2y = 0, l2z = 0.2025;
    l3x = 0, l3y = 0.2045, l3z = 0;
    l4x = 0, l4y = 0, l4z = 0.2155;
    l5x = 0, l5y = 0.1845, l5z = 0;
    l6x = 0, l6y = 0, l6z = 0.2155;
    l7x = 0,l7y = 0.081, l7z = 0;
    
    m1 = 5.76; m2 = 6.35; m3 = 3.5; m4 = 3.5; m5 = 3.5; m6 = 1.8; m7 = 1.2;
    g = 9.81;
    
    % viscous coefficients of positive and negative directions
    bv1_positive=0.1392;  bv2_positive=8.3008; bv3_positive=0.0108;  bv4_positive=1.4523; bv5_positive=0.0178;  bv6_positive=0.1199; bv7_positive= -0.0573;
    bv1_negative=0.2439;  bv2_negative=8.6682; bv3_negative=-0.0106;  bv4_negative=1.6164; bv5_negative=0.0104;  bv6_negative=0.0961; bv7_negative= -0.0319;
    % coulomb coefficients of positive and negative directions
    bc1_positive=1.7809;  bc2_positive= 1.9103; bc3_positive= -0.0445;  bc4_positive=0.2538; bc5_positive=0.1151;  bc6_positive=0.0534; bc7_positive=0.4934;
    bc1_negative= -0.3484;  bc2_negative= 0.1029; bc3_negative= -0.3254;  bc4_negative=-0.4345; bc5_negative=0.0809;  bc6_negative=-0.0881; bc7_negative=-0.1781;
    %    bv1_positive=0; bv2_positive=0;
    %lc1 = .5; lc2 = 1;
    
    c1x = 0, c1y = -0.03, c1z = 0.12;
    c2x = 0.0003, c2y = 0.059, c2z = 0.042;
    c3x = 0, c3y = 0.03, c3z = 0.13;
    c4x = 0, c4y = 0.067, c4z = 0.034;
    c5x = 0.0001, c5y = 0.021, c5z = 0.076;
    c6x = 0, c6y = 0.0006, c6z = 0.0004;
    c7x = 0, c7y = 0, c7z = 0.02;
    
    I1xx= 0.033, I1xy= 0, I1xz= 0, I1yy= 0.0333, I1yz= 0.004887, I1zz= 0.0123;
    I2xx= 0.0305, I2xy= 0, I2xz= 0, I2yy= 0.0304, I2yz= 0.004887, I2zz= 0.011;
    I3xx= 0.025, I3xy= 0, I3xz= 0, I3yy= 0.0238, I3yz= 0.00487, I3zz= 0.0076;
    I4xx= 0.017, I4xy= 0, I4xz= 0, I4yy= 0.0164, I4yz= 0.00284, I4zz= 0.006;
    I5xx= 0.01, I5xy= 0, I5xz= 0, I5yy= 0.0087, I5yz= 0.00309, I5zz= 0.00449;
    I6xx= 0.0049, I6xy= 0, I6xz= 0, I6yy= 0.0047, I6yz= 0.000246, I6zz= 0.0036;
    I7xx= 0.0002, I7xy= 0, I7xz= 0, I7yy= 0.0002, I7yz= 0, I7zz= 0.0003;
    
    xG
    uG
  end
  
  methods
    function obj = KukaArmPlant
      obj = obj@Manipulator(7,7);%[To Be Checked]
      obj = setInputLimits(obj,-10,10);%[To Be Checked]

      obj = setInputFrame(obj,CoordinateFrame('KukaArmInput',7,'u',{'tau1','tau2','tau3','tau4','tau5','tau6','tau7'}));
      obj = setStateFrame(obj,CoordinateFrame('KukaArmState',14,'x',{'theta1','theta2','theta3','theta4','theta5','theta6','theta7', ...
                                                'theta1dot','theta2dot','theta3dot','theta4dot','theta5dot','theta6dot','theta7dot'}));
      obj = setOutputFrame(obj,obj.getStateFrame);
      
      obj.xG = Point(obj.getStateFrame,[pi;pi/2;pi/2;pi/2;pi/2;pi/2;pi/2;0;0;0;0;0;0;0]);%[To Be Checked]
      obj.uG = Point(obj.getInputFrame,[0;0;0;0;0;0;0]);

%       obj = setParamFrame(obj,CoordinateFrame('KukaArmParams',28,'p',{'m1','m2','m3','m4','m5','m6','m7','I1xx','I1yy','I1zz','I2xx', ...
%           'I2yy','I2zz','I3xx','I3yy','I3zz','I4xx','I4yy','I4zz','I5xx','I5yy','I5zz','I6xx','I6yy','I6zz','I7xx','I7yy','I7zz'}));
%       obj = setParamFrame(obj,CoordinateFrame('KukaArmParams',14,'p',...
%         {'bv1_positive','bv2_positive','bv3_positive','bv4_positive','bv5_positive','bv6_positive','bv7_positive','m1','m2','m3','m4','m5','m6','m7'}));
      obj = setParamFrame(obj,CoordinateFrame('KukaArmParams',1,'p',...
        {'m1'})); %,'m2','m3','m4','m5','m6','m7'
%       obj = setParamFrame(obj,CoordinateFrame('KukaArmParams',6,'p',...
%         { 'bv1_positive','bv2_positive','lc1','lc2','Ic1','Ic2' }));
%      obj = setParamFrame(obj,CoordinateFrame('KukaArmParams',10,'p',...
%        { 'l1','l2','m1','m2','bv1_positive','bv2_positive','lc1','lc2','I1','I2' }));
%      obj = setParamFrame(obj,CoordinateFrame('KukaArmParams',1,'p',...
%        { 'lc2' }));
      obj = setParamLimits(obj,zeros(obj.getParamFrame.dim,1));
    end

    function [H,C,B] = manipulatorDynamics(obj,q,qd)
      % keep it readable:
      m1=obj.m1; m2=obj.m2; m3=obj.m3; m4=obj.m4; m5=obj.m5; m6=obj.m6; m7=obj.m7;
      l1x=obj.l1x; l1y=obj.l1y; l1z=obj.l1z;
      l2x=obj.l2x; l2y=obj.l2y; l2z=obj.l2z;
      l3x=obj.l3x; l3y=obj.l3y; l3z=obj.l3z;
      l4x=obj.l4x; l4y=obj.l4y; l4z=obj.l4z;
      l5x=obj.l5x; l5y=obj.l5y; l5z=obj.l5z;
      l6x=obj.l6x; l6y=obj.l6y; l6z=obj.l6z;
      l7x=obj.l7x; l7y=obj.l7y; l7z=obj.l7z;
      
      g=obj.g; 
      c1x = obj.c1x; c1y = obj.c1y; c1z = obj.c1z;
      c2x = obj.c2x; c2y = obj.c2y; c2z = obj.c2z;
      c3x = obj.c3x; c3y = obj.c3y; c3z = obj.c3z;
      c4x = obj.c4x; c4y = obj.c4y; c4z = obj.c4z;
      c5x = obj.c5x; c5y = obj.c5y; c5z = obj.c5z;
      c6x = obj.c6x; c6y = obj.c6y; c6z = obj.c6z;
      c7x = obj.c7x; c7y = obj.c7y; c7z = obj.c7z;
      
      bv1_positive=obj.bv1_positive; bv2_positive=obj.bv2_positive; bv3_positive=obj.bv3_positive; bv4_positive=obj.bv4_positive; bv5_positive=obj.bv5_positive; bv6_positive=obj.bv6_positive; bv7_positive=obj.bv7_positive;
      bv1_negative=obj.bv1_negative; bv2_negative=obj.bv2_negative; bv3_negative=obj.bv3_negative; bv4_negative=obj.bv4_negative; bv5_negative=obj.bv5_negative; bv6_negative=obj.bv6_negative; bv7_negative=obj.bv7_negative;
      bc1_positive=obj.bc1_positive; bc2_positive=obj.bc2_positive; bc3_positive=obj.bc3_positive; bc4_positive=obj.bc4_positive; bc5_positive=obj.bc5_positive; bc6_positive=obj.bc6_positive; bc7_positive=obj.bc7_positive;
      bc1_negative=obj.bc1_negative; bc2_negative=obj.bc2_negative; bc3_negative=obj.bc3_negative; bc4_negative=obj.bc4_negative; bc5_negative=obj.bc5_negative; bc6_negative=obj.bc6_negative; bc7_negative=obj.bc7_negative;
      
      %I1 = obj.Ic1 + obj.m1*obj.lc1^2; I2 = obj.Ic2 + obj.m2*obj.lc2^2;
      %I1xx = obj.I1xx; I2xx = obj.I2xx;
      %m2l1lc2 = m2*l1*lc2;  % occurs often!
      tic
      c = cos(q(1:7,:));  s = sin(q(1:7,:));%[To be Checked]
      toc
      
      %h12 = I2 + m2l1lc2*c(2);
      %h11 = I1zz + (sin(q2)*(sin(q3)*sin(q5) + cos(q3)*cos(q4)*cos(q5)) + cos(q2)*cos(q5)*sin(q4))*(I5xz*(cos(q2)*cos(q4) - cos(q3)*sin(q2)*sin(q4)) + I5xx*(sin(q2)*(sin(q3)*sin(q5) + cos(q3)*cos(q4)*cos(q5)) + cos(q2)*cos(q5)*sin(q4)) + I5xy*(sin(q2)*(cos(q5)*sin(q3) - cos(q3)*cos(q4)*sin(q5)) - cos(q2)*sin(q4)*sin(q5))) + (cos(q2)*(sin(q4)*(sin(q5)*sin(q7) + cos(q5)*cos(q6)*cos(q7)) + cos(q4)*cos(q7)*sin(q6)) + sin(q2)*(cos(q3)*(cos(q4)*(sin(q5)*sin(q7) + cos(q5)*cos(q6)*cos(q7)) - cos(q7)*sin(q4)*sin(q6)) - sin(q3)*(cos(q5)*sin(q7) - cos(q6)*cos(q7)*sin(q5))))*(I7xx*(cos(q2)*(sin(q4)*(sin(q5)*sin(q7) + cos(q5)*cos(q6)*cos(q7)) + cos(q4)*cos(q7)*sin(q6)) + sin(q2)*(cos(q3)*(cos(q4)*(sin(q5)*sin(q7) + cos(q5)*cos(q6)*cos(q7)) - cos(q7)*sin(q4)*sin(q6)) - sin(q3)*(cos(q5)*sin(q7) - cos(q6)*cos(q7)*sin(q5)))) + I7xy*(cos(q2)*(sin(q4)*(cos(q7)*sin(q5) - cos(q5)*cos(q6)*sin(q7)) - cos(q4)*sin(q6)*sin(q7)) + sin(q2)*(cos(q3)*(cos(q4)*(cos(q7)*sin(q5) - cos(q5)*cos(q6)*sin(q7)) + sin(q4)*sin(q6)*sin(q7)) - sin(q3)*(cos(q5)*cos(q7) + cos(q6)*sin(q5)*sin(q7)))) - I7xz*(sin(q2)*(cos(q3)*(cos(q6)*sin(q4) + cos(q4)*cos(q5)*sin(q6)) + sin(q3)*sin(q5)*sin(q6)) - cos(q2)*(cos(q4)*cos(q6) - cos(q5)*sin(q4)*sin(q6)))) + (sin(q2)*(cos(q5)*sin(q3) - cos(q3)*cos(q4)*sin(q5)) - cos(q2)*sin(q4)*sin(q5))*(I5yz*(cos(q2)*cos(q4) - cos(q3)*sin(q2)*sin(q4)) + I5xy*(sin(q2)*(sin(q3)*sin(q5) + cos(q3)*cos(q4)*cos(q5)) + cos(q2)*cos(q5)*sin(q4)) + I5yy*(sin(q2)*(cos(q5)*sin(q3) - cos(q3)*cos(q4)*sin(q5)) - cos(q2)*sin(q4)*sin(q5))) + m5*((l1y*cos(q1) + l1x*sin(q1) - sin(q1)*(l2x*cos(q2) - l2y*sin(q2) + sin(q2)*(l3z - c5z*cos(q4) + l4y*cos(q4) + l4x*sin(q4) - sin(q4)*(c5x*cos(q5) - c5y*sin(q5))) - cos(q2)*(cos(q3)*(l4x*cos(q4) + c5z*sin(q4) - l4y*sin(q4) - cos(q4)*(c5x*cos(q5) - c5y*sin(q5))) + l3x*cos(q3) - sin(q3)*(c5y*cos(q5) - l4z + c5x*sin(q5)) - l3y*sin(q3))) + cos(q1)*(sin(q3)*(l4x*cos(q4) + c5z*sin(q4) - l4y*sin(q4) - cos(q4)*(c5x*cos(q5) - c5y*sin(q5))) - l2z + cos(q3)*(c5y*cos(q5) - l4z + c5x*sin(q5)) + l3y*cos(q3) + l3x*sin(q3)))^2 + (cos(q1)*(l2x*cos(q2) - l2y*sin(q2) + sin(q2)*(l3z - c5z*cos(q4) + l4y*cos(q4) + l4x*sin(q4) - sin(q4)*(c5x*cos(q5) - c5y*sin(q5))) - cos(q2)*(cos(q3)*(l4x*cos(q4) + c5z*sin(q4) - l4y*sin(q4) - cos(q4)*(c5x*cos(q5) - c5y*sin(q5))) + l3x*cos(q3) - sin(q3)*(c5y*cos(q5) - l4z + c5x*sin(q5)) - l3y*sin(q3))) - l1x*cos(q1) + l1y*sin(q1) + sin(q1)*(sin(q3)*(l4x*cos(q4) + c5z*sin(q4) - l4y*sin(q4) - cos(q4)*(c5x*cos(q5) - c5y*sin(q5))) - l2z + cos(q3)*(c5y*cos(q5) - l4z + c5x*sin(q5)) + l3y*cos(q3) + l3x*sin(q3)))^2) + m1*((c1y*cos(q1) + c1x*sin(q1))^2 + (c1x*cos(q1) - c1y*sin(q1))^2) + (cos(q2)*(sin(q4)*(cos(q7)*sin(q5) - cos(q5)*cos(q6)*sin(q7)) - cos(q4)*sin(q6)*sin(q7)) + sin(q2)*(cos(q3)*(cos(q4)*(cos(q7)*sin(q5) - cos(q5)*cos(q6)*sin(q7)) + sin(q4)*sin(q6)*sin(q7)) - sin(q3)*(cos(q5)*cos(q7) + cos(q6)*sin(q5)*sin(q7))))*(I7xy*(cos(q2)*(sin(q4)*(sin(q5)*sin(q7) + cos(q5)*cos(q6)*cos(q7)) + cos(q4)*cos(q7)*sin(q6)) + sin(q2)*(cos(q3)*(cos(q4)*(sin(q5)*sin(q7) + cos(q5)*cos(q6)*cos(q7)) - cos(q7)*sin(q4)*sin(q6)) - sin(q3)*(cos(q5)*sin(q7) - cos(q6)*cos(q7)*sin(q5)))) + I7yy*(cos(q2)*(sin(q4)*(cos(q7)*sin(q5) - cos(q5)*cos(q6)*sin(q7)) - cos(q4)*sin(q6)*sin(q7)) + sin(q2)*(cos(q3)*(cos(q4)*(cos(q7)*sin(q5) - cos(q5)*cos(q6)*sin(q7)) + sin(q4)*sin(q6)*sin(q7)) - sin(q3)*(cos(q5)*cos(q7) + cos(q6)*sin(q5)*sin(q7)))) - I7yz*(sin(q2)*(cos(q3)*(cos(q6)*sin(q4) + cos(q4)*cos(q5)*sin(q6)) + sin(q3)*sin(q5)*sin(q6)) - cos(q2)*(cos(q4)*cos(q6) - cos(q5)*sin(q4)*sin(q6)))) + (sin(q2)*(cos(q5)*sin(q3) - cos(q3)*cos(q4)*sin(q5)) - cos(q2)*sin(q4)*sin(q5))*(I6zz*(sin(q2)*(cos(q5)*sin(q3) - cos(q3)*cos(q4)*sin(q5)) - cos(q2)*sin(q4)*sin(q5)) + I6xz*(sin(q2)*(cos(q3)*(sin(q4)*sin(q6) - cos(q4)*cos(q5)*cos(q6)) - cos(q6)*sin(q3)*sin(q5)) - cos(q2)*(cos(q4)*sin(q6) + cos(q5)*cos(q6)*sin(q4))) + I6yz*(sin(q2)*(cos(q3)*(cos(q6)*sin(q4) + cos(q4)*cos(q5)*sin(q6)) + sin(q3)*sin(q5)*sin(q6)) - cos(q2)*(cos(q4)*cos(q6) - cos(q5)*sin(q4)*sin(q6)))) + m2*((c2z*cos(q1) - l1y*cos(q1) - l1x*sin(q1) + sin(q1)*(c2x*cos(q2) - c2y*sin(q2)))^2 + (l1x*cos(q1) + c2z*sin(q1) - l1y*sin(q1) - cos(q1)*(c2x*cos(q2) - c2y*sin(q2)))^2) + m7*((cos(q1)*(l3y*cos(q3) - l2z + l3x*sin(q3) + sin(q3)*(l4x*cos(q4) - l4y*sin(q4) + sin(q4)*(l5z - c7z*cos(q6) + l6y*cos(q6) + l6x*sin(q6) - sin(q6)*(c7x*cos(q7) - c7y*sin(q7))) - cos(q4)*(cos(q5)*(l6x*cos(q6) + c7z*sin(q6) - l6y*sin(q6) - cos(q6)*(c7x*cos(q7) - c7y*sin(q7))) + l5x*cos(q5) - sin(q5)*(c7y*cos(q7) - l6z + c7x*sin(q7)) - l5y*sin(q5))) + cos(q3)*(sin(q5)*(l6x*cos(q6) + c7z*sin(q6) - l6y*sin(q6) - cos(q6)*(c7x*cos(q7) - c7y*sin(q7))) - l4z + cos(q5)*(c7y*cos(q7) - l6z + c7x*sin(q7)) + l5y*cos(q5) + l5x*sin(q5))) - sin(q1)*(sin(q2)*(l3z + l4y*cos(q4) + l4x*sin(q4) - cos(q4)*(l5z - c7z*cos(q6) + l6y*cos(q6) + l6x*sin(q6) - sin(q6)*(c7x*cos(q7) - c7y*sin(q7))) - sin(q4)*(cos(q5)*(l6x*cos(q6) + c7z*sin(q6) - l6y*sin(q6) - cos(q6)*(c7x*cos(q7) - c7y*sin(q7))) + l5x*cos(q5) - sin(q5)*(c7y*cos(q7) - l6z + c7x*sin(q7)) - l5y*sin(q5))) - cos(q2)*(l3x*cos(q3) + cos(q3)*(l4x*cos(q4) - l4y*sin(q4) + sin(q4)*(l5z - c7z*cos(q6) + l6y*cos(q6) + l6x*sin(q6) - sin(q6)*(c7x*cos(q7) - c7y*sin(q7))) - cos(q4)*(cos(q5)*(l6x*cos(q6) + c7z*sin(q6) - l6y*sin(q6) - cos(q6)*(c7x*cos(q7) - c7y*sin(q7))) + l5x*cos(q5) - sin(q5)*(c7y*cos(q7) - l6z + c7x*sin(q7)) - l5y*sin(q5))) - l3y*sin(q3) - sin(q3)*(sin(q5)*(l6x*cos(q6) + c7z*sin(q6) - l6y*sin(q6) - cos(q6)*(c7x*cos(q7) - c7y*sin(q7))) - l4z + cos(q5)*(c7y*cos(q7) - l6z + c7x*sin(q7)) + l5y*cos(q5) + l5x*sin(q5))) + l2x*cos(q2) - l2y*sin(q2)) + l1y*cos(q1) + l1x*sin(q1))^2 + (sin(q1)*(l3y*cos(q3) - l2z + l3x*sin(q3) + sin(q3)*(l4x*cos(q4) - l4y*sin(q4) + sin(q4)*(l5z - c7z*cos(q6) + l6y*cos(q6) + l6x*sin(q6) - sin(q6)*(c7x*cos(q7) - c7y*sin(q7))) - cos(q4)*(cos(q5)*(l6x*cos(q6) + c7z*sin(q6) - l6y*sin(q6) - cos(q6)*(c7x*cos(q7) - c7y*sin(q7))) + l5x*cos(q5) - sin(q5)*(c7y*cos(q7) - l6z + c7x*sin(q7)) - l5y*sin(q5))) + cos(q3)*(sin(q5)*(l6x*cos(q6) + c7z*sin(q6) - l6y*sin(q6) - cos(q6)*(c7x*cos(q7) - c7y*sin(q7))) - l4z + cos(q5)*(c7y*cos(q7) - l6z + c7x*sin(q7)) + l5y*cos(q5) + l5x*sin(q5))) + cos(q1)*(sin(q2)*(l3z + l4y*cos(q4) + l4x*sin(q4) - cos(q4)*(l5z - c7z*cos(q6) + l6y*cos(q6) + l6x*sin(q6) - sin(q6)*(c7x*cos(q7) - c7y*sin(q7))) - sin(q4)*(cos(q5)*(l6x*cos(q6) + c7z*sin(q6) - l6y*sin(q6) - cos(q6)*(c7x*cos(q7) - c7y*sin(q7))) + l5x*cos(q5) - sin(q5)*(c7y*cos(q7) - l6z + c7x*sin(q7)) - l5y*sin(q5))) - cos(q2)*(l3x*cos(q3) + cos(q3)*(l4x*cos(q4) - l4y*sin(q4) + sin(q4)*(l5z - c7z*cos(q6) + l6y*cos(q6) + l6x*sin(q6) - sin(q6)*(c7x*cos(q7) - c7y*sin(q7))) - cos(q4)*(cos(q5)*(l6x*cos(q6) + c7z*sin(q6) - l6y*sin(q6) - cos(q6)*(c7x*cos(q7) - c7y*sin(q7))) + l5x*cos(q5) - sin(q5)*(c7y*cos(q7) - l6z + c7x*sin(q7)) - l5y*sin(q5))) - l3y*sin(q3) - sin(q3)*(sin(q5)*(l6x*cos(q6) + c7z*sin(q6) - l6y*sin(q6) - cos(q6)*(c7x*cos(q7) - c7y*sin(q7))) - l4z + cos(q5)*(c7y*cos(q7) - l6z + c7x*sin(q7)) + l5y*cos(q5) + l5x*sin(q5))) + l2x*cos(q2) - l2y*sin(q2)) - l1x*cos(q1) + l1y*sin(q1))^2) - (sin(q2)*(cos(q3)*(cos(q6)*sin(q4) + cos(q4)*cos(q5)*sin(q6)) + sin(q3)*sin(q5)*sin(q6)) - cos(q2)*(cos(q4)*cos(q6) - cos(q5)*sin(q4)*sin(q6)))*(I7xz*(cos(q2)*(sin(q4)*(sin(q5)*sin(q7) + cos(q5)*cos(q6)*cos(q7)) + cos(q4)*cos(q7)*sin(q6)) + sin(q2)*(cos(q3)*(cos(q4)*(sin(q5)*sin(q7) + cos(q5)*cos(q6)*cos(q7)) - cos(q7)*sin(q4)*sin(q6)) - sin(q3)*(cos(q5)*sin(q7) - cos(q6)*cos(q7)*sin(q5)))) + I7yz*(cos(q2)*(sin(q4)*(cos(q7)*sin(q5) - cos(q5)*cos(q6)*sin(q7)) - cos(q4)*sin(q6)*sin(q7)) + sin(q2)*(cos(q3)*(cos(q4)*(cos(q7)*sin(q5) - cos(q5)*cos(q6)*sin(q7)) + sin(q4)*sin(q6)*sin(q7)) - sin(q3)*(cos(q5)*cos(q7) + cos(q6)*sin(q5)*sin(q7)))) - I7zz*(sin(q2)*(cos(q3)*(cos(q6)*sin(q4) + cos(q4)*cos(q5)*sin(q6)) + sin(q3)*sin(q5)*sin(q6)) - cos(q2)*(cos(q4)*cos(q6) - cos(q5)*sin(q4)*sin(q6)))) + m3*((cos(q1)*(l2x*cos(q2) + c3z*sin(q2) - l2y*sin(q2) - cos(q2)*(c3x*cos(q3) - c3y*sin(q3))) - l1x*cos(q1) + sin(q1)*(c3y*cos(q3) - l2z + c3x*sin(q3)) + l1y*sin(q1))^2 + (cos(q1)*(c3y*cos(q3) - l2z + c3x*sin(q3)) - sin(q1)*(l2x*cos(q2) + c3z*sin(q2) - l2y*sin(q2) - cos(q2)*(c3x*cos(q3) - c3y*sin(q3))) + l1y*cos(q1) + l1x*sin(q1))^2) + (sin(q2)*(cos(q3)*(sin(q4)*sin(q6) - cos(q4)*cos(q5)*cos(q6)) - cos(q6)*sin(q3)*sin(q5)) - cos(q2)*(cos(q4)*sin(q6) + cos(q5)*cos(q6)*sin(q4)))*(I6xz*(sin(q2)*(cos(q5)*sin(q3) - cos(q3)*cos(q4)*sin(q5)) - cos(q2)*sin(q4)*sin(q5)) + I6xx*(sin(q2)*(cos(q3)*(sin(q4)*sin(q6) - cos(q4)*cos(q5)*cos(q6)) - cos(q6)*sin(q3)*sin(q5)) - cos(q2)*(cos(q4)*sin(q6) + cos(q5)*cos(q6)*sin(q4))) + I6xy*(sin(q2)*(cos(q3)*(cos(q6)*sin(q4) + cos(q4)*cos(q5)*sin(q6)) + sin(q3)*sin(q5)*sin(q6)) - cos(q2)*(cos(q4)*cos(q6) - cos(q5)*sin(q4)*sin(q6)))) + m6*((cos(q1)*(cos(q3)*(l5y*cos(q5) - c6z*cos(q5) - l4z + l5x*sin(q5) + sin(q5)*(c6x*cos(q6) - c6y*sin(q6))) - l2z + sin(q3)*(sin(q4)*(l5z + c6y*cos(q6) + c6x*sin(q6)) - cos(q4)*(l5x*cos(q5) + c6z*sin(q5) - l5y*sin(q5) + cos(q5)*(c6x*cos(q6) - c6y*sin(q6))) + l4x*cos(q4) - l4y*sin(q4)) + l3y*cos(q3) + l3x*sin(q3)) - sin(q1)*(sin(q2)*(l3z - cos(q4)*(l5z + c6y*cos(q6) + c6x*sin(q6)) - sin(q4)*(l5x*cos(q5) + c6z*sin(q5) - l5y*sin(q5) + cos(q5)*(c6x*cos(q6) - c6y*sin(q6))) + l4y*cos(q4) + l4x*sin(q4)) + l2x*cos(q2) + cos(q2)*(sin(q3)*(l5y*cos(q5) - c6z*cos(q5) - l4z + l5x*sin(q5) + sin(q5)*(c6x*cos(q6) - c6y*sin(q6))) - cos(q3)*(sin(q4)*(l5z + c6y*cos(q6) + c6x*sin(q6)) - cos(q4)*(l5x*cos(q5) + c6z*sin(q5) - l5y*sin(q5) + cos(q5)*(c6x*cos(q6) - c6y*sin(q6))) + l4x*cos(q4) - l4y*sin(q4)) - l3x*cos(q3) + l3y*sin(q3)) - l2y*sin(q2)) + l1y*cos(q1) + l1x*sin(q1))^2 + (sin(q1)*(cos(q3)*(l5y*cos(q5) - c6z*cos(q5) - l4z + l5x*sin(q5) + sin(q5)*(c6x*cos(q6) - c6y*sin(q6))) - l2z + sin(q3)*(sin(q4)*(l5z + c6y*cos(q6) + c6x*sin(q6)) - cos(q4)*(l5x*cos(q5) + c6z*sin(q5) - l5y*sin(q5) + cos(q5)*(c6x*cos(q6) - c6y*sin(q6))) + l4x*cos(q4) - l4y*sin(q4)) + l3y*cos(q3) + l3x*sin(q3)) + cos(q1)*(sin(q2)*(l3z - cos(q4)*(l5z + c6y*cos(q6) + c6x*sin(q6)) - sin(q4)*(l5x*cos(q5) + c6z*sin(q5) - l5y*sin(q5) + cos(q5)*(c6x*cos(q6) - c6y*sin(q6))) + l4y*cos(q4) + l4x*sin(q4)) + l2x*cos(q2) + cos(q2)*(sin(q3)*(l5y*cos(q5) - c6z*cos(q5) - l4z + l5x*sin(q5) + sin(q5)*(c6x*cos(q6) - c6y*sin(q6))) - cos(q3)*(sin(q4)*(l5z + c6y*cos(q6) + c6x*sin(q6)) - cos(q4)*(l5x*cos(q5) + c6z*sin(q5) - l5y*sin(q5) + cos(q5)*(c6x*cos(q6) - c6y*sin(q6))) + l4x*cos(q4) - l4y*sin(q4)) - l3x*cos(q3) + l3y*sin(q3)) - l2y*sin(q2)) - l1x*cos(q1) + l1y*sin(q1))^2) + (sin(q2)*(cos(q3)*(cos(q6)*sin(q4) + cos(q4)*cos(q5)*sin(q6)) + sin(q3)*sin(q5)*sin(q6)) - cos(q2)*(cos(q4)*cos(q6) - cos(q5)*sin(q4)*sin(q6)))*(I6yz*(sin(q2)*(cos(q5)*sin(q3) - cos(q3)*cos(q4)*sin(q5)) - cos(q2)*sin(q4)*sin(q5)) + I6xy*(sin(q2)*(cos(q3)*(sin(q4)*sin(q6) - cos(q4)*cos(q5)*cos(q6)) - cos(q6)*sin(q3)*sin(q5)) - cos(q2)*(cos(q4)*sin(q6) + cos(q5)*cos(q6)*sin(q4))) + I6yy*(sin(q2)*(cos(q3)*(cos(q6)*sin(q4) + cos(q4)*cos(q5)*sin(q6)) + sin(q3)*sin(q5)*sin(q6)) - cos(q2)*(cos(q4)*cos(q6) - cos(q5)*sin(q4)*sin(q6)))) + (cos(q2)*cos(q4) - cos(q3)*sin(q2)*sin(q4))*(I5zz*(cos(q2)*cos(q4) - cos(q3)*sin(q2)*sin(q4)) + I5xz*(sin(q2)*(sin(q3)*sin(q5) + cos(q3)*cos(q4)*cos(q5)) + cos(q2)*cos(q5)*sin(q4)) + I5yz*(sin(q2)*(cos(q5)*sin(q3) - cos(q3)*cos(q4)*sin(q5)) - cos(q2)*sin(q4)*sin(q5))) + cos(q2)*(I2yy*cos(q2) + I2xy*sin(q2)) + sin(q2)*(I2xy*cos(q2) + I2xx*sin(q2)) + (cos(q2)*sin(q4) + cos(q3)*cos(q4)*sin(q2))*(I4xx*(cos(q2)*sin(q4) + cos(q3)*cos(q4)*sin(q2)) + I4xy*(cos(q2)*cos(q4) - cos(q3)*sin(q2)*sin(q4)) + I4xz*sin(q2)*sin(q3)) + (cos(q2)*cos(q4) - cos(q3)*sin(q2)*sin(q4))*(I4xy*(cos(q2)*sin(q4) + cos(q3)*cos(q4)*sin(q2)) + I4yy*(cos(q2)*cos(q4) - cos(q3)*sin(q2)*sin(q4)) + I4yz*sin(q2)*sin(q3)) + m4*((cos(q1)*(l3y*cos(q3) - c4z*cos(q3) - l2z + l3x*sin(q3) + sin(q3)*(c4x*cos(q4) - c4y*sin(q4))) - sin(q1)*(sin(q2)*(l3z + c4y*cos(q4) + c4x*sin(q4)) - cos(q2)*(l3x*cos(q3) + c4z*sin(q3) - l3y*sin(q3) + cos(q3)*(c4x*cos(q4) - c4y*sin(q4))) + l2x*cos(q2) - l2y*sin(q2)) + l1y*cos(q1) + l1x*sin(q1))^2 + (sin(q1)*(l3y*cos(q3) - c4z*cos(q3) - l2z + l3x*sin(q3) + sin(q3)*(c4x*cos(q4) - c4y*sin(q4))) + cos(q1)*(sin(q2)*(l3z + c4y*cos(q4) + c4x*sin(q4)) - cos(q2)*(l3x*cos(q3) + c4z*sin(q3) - l3y*sin(q3) + cos(q3)*(c4x*cos(q4) - c4y*sin(q4))) + l2x*cos(q2) - l2y*sin(q2)) - l1x*cos(q1) + l1y*sin(q1))^2) + cos(q2)*(I3zz*cos(q2) + I3xz*cos(q3)*sin(q2) - I3yz*sin(q2)*sin(q3)) + sin(q2)*sin(q3)*(I4xz*(cos(q2)*sin(q4) + cos(q3)*cos(q4)*sin(q2)) + I4yz*(cos(q2)*cos(q4) - cos(q3)*sin(q2)*sin(q4)) + I4zz*sin(q2)*sin(q3)) + cos(q3)*sin(q2)*(I3xz*cos(q2) + I3xx*cos(q3)*sin(q2) - I3xy*sin(q2)*sin(q3)) - sin(q2)*sin(q3)*(I3yz*cos(q2) + I3xy*cos(q3)*sin(q2) - I3yy*sin(q2)*sin(q3));
      %H = [ I1 + I2 + m2*l1^2 + 2*m2l1lc2*c(2), h12; h12, I2 ];
      
      tic
      H = getKukaArmInertiaMatrix_simplified(obj,q,c,s);
      toc
      
      tic
      C = getKukaArmCoriolisCentrifugalVector(obj,q,qd,c,s);
      toc
      %C = [ -2*m2l1lc2*s(2)*qd(2), -m2l1lc2*s(2)*qd(2); m2l1lc2*s(2)*qd(1), 0 ];
      %G = g*[ m1*lc1*s(1) + m2*(l1*s(1)+lc2*s12); m2*lc2*s12 ];
      tic
      G = getKukaArmGravityVector(obj,q,c,s);
      toc
      % Coulomb and viscous friction
      
%       for i =1:7
%         if (qd(i) > 0)
%             b = bv1_positive * qd(i) + bc1_positive; 
%         else
%             b = bv1_negative * qd(i) + bc1_negative;
%         end
%       end

      for i =1:7
          b(i,1) = (bv1_positive + bv1_positive)/2 * qd(i) + bc1_positive;
      end
      
      % accumate total C and add a damping term:
      C = C + G + b;
      B = eye(7);
    end
    
    function [T,U] = energy(obj,x)
      m1=obj.m1; m2=obj.m2; m3=obj.m3; m4=obj.m4; m5=obj.m5; m6=obj.m6; m7=obj.m7; 
      
      l1x=obj.l1x; l1y=obj.l1y; l1z=obj.l1z;
      l2x=obj.l2x; l2y=obj.l2y; l2z=obj.l2z;
      l3x=obj.l3x; l3y=obj.l3y; l3z=obj.l3z;
      l4x=obj.l4x; l4y=obj.l4y; l4z=obj.l4z;
      l5x=obj.l5x; l5y=obj.l5y; l5z=obj.l5z;
      l6x=obj.l6x; l6y=obj.l6y; l6z=obj.l6z;
      l7x=obj.l7x; l7y=obj.l7y; l7z=obj.l7z;
      
      g=obj.g; 
      c1x = obj.c1x; c1y = obj.c1y; c1z = obj.c1z;
      c2x = obj.c2x; c2y = obj.c2y; c2z = obj.c2z;
      c3x = obj.c3x; c3y = obj.c3y; c3z = obj.c3z;
      c4x = obj.c4x; c4y = obj.c4y; c4z = obj.c4z;
      c5x = obj.c5x; c5y = obj.c5y; c5z = obj.c5z;
      c6x = obj.c6x; c6y = obj.c6y; c6z = obj.c6z;
      c7x = obj.c7x; c7y = obj.c7y; c7z = obj.c7z;
      
%       bv1_positive=obj.bv1_positive; bv2_positive=obj.bv2_positive; bv3_positive=obj.bv3_positive; bv4_positive=obj.bv4_positive; bv5_positive=obj.bv5_positive; bv6_positive=obj.bv6_positive; bv7_positive=obj.bv7_positive;
%       bv1_negative=obj.bv1_negative; bv2_negative=obj.bv2_negative; bv3_negative=obj.bv3_negative; bv4_negative=obj.bv4_negative; bv5_negative=obj.bv5_negative; bv6_negative=obj.bv6_negative; bv7_negative=obj.bv7_negative;
%       bc1_positive=obj.bc1_positive; bc2_positive=obj.bc2_positive; bc3_positive=obj.bc3_positive; bc4_positive=obj.bc4_positive; bc5_positive=obj.bc5_positive; bc6_positive=obj.bc6_positive; bc7_positive=obj.bc7_positive;
%       bc1_negative=obj.bc1_negative; bc2_negative=obj.bc2_negative; bc3_negative=obj.bc3_negative; bc4_negative=obj.bc4_negative; bc5_negative=obj.bc5_negative; bc6_negative=obj.bc6_negative; bc7_negative=obj.bc7_negative;
      
      %I1 = obj.Ic1 + obj.m1*obj.lc1^2; I2 = obj.Ic2 + obj.m2*obj.lc2^2;
      %I1xx = obj.I1xx; I2xx = obj.I2xx;
      q = x(1:7); qd = x(8:14);
      c = cos(q(1:7,:));  s = sin(q(1:7,:));
      
%       T = .5*I1*qd(1)^2 + .5*(m2*l1^2 + I2 + 2*m2*l1*lc2*c(2))*qd(1)^2 + .5*I2*qd(2)^2 + (I2 + m2*l1*lc2*c(2))*qd(1)*qd(2);
%       U = -m1*g*lc1*c(1) - m2*g*(l1*c(1)+lc2*c12);      
      %T = .5*I1xx*qd(1)^2 + .5*(m2*l1x^2 + I2xx + 2*m2*lx1*c2x*c(2))*qd(1)^2 + .5*I2xx*qd(2)^2 + (I2xx + m2*l1x*c2x*c(2))*qd(1)*qd(2);
      %KineticSymbol = fopen('KineticEnergy.txt');
      %T = textscan(KineticSymbol);
      %fclose(KineticSymbol);
      %T = readfile('KineticEnergy.txt'); %getKineticEnergy();
      
      T = getKineticEnergy(obj,q,qd);
      U = getPotentialEnergy(obj,q);
      
      %U = -m1*g*c1x*c(1) - m2*g*(l1x*c(1)+c2x*c12);
    end
    
    % todo: also implement sodynamics here so that I can keep the
    % vectorized version?

%[To be checked]
%     function [f,df,d2f,d3f] = dynamics(obj,t,x,u)
%       f = dynamics@Manipulator(obj,t,x,u);
%       if (nargout>1)
%         [df,d2f,d3f]= dynamicsGradients(obj,t,x,u,nargout-1);
%       end
%     end
    
    function x = getInitialState(obj)
      x = .1*randn(14,1);
    end
    
    function n = getNumPositions(obj)
      n = 7;
    end
    
    function n = getNumVelocities(obj)
      n = 7;
    end
    
    function [c,V]=balanceLQR(obj)
      Q = diag([10,10,1,1]); R = 1;
      if (nargout<2)
        c = tilqr(obj,obj.xG,obj.uG,Q,R);
      else
        if any(~isinf([obj.umin;obj.umax]))
          error('currently, you must disable input limits to estimate the ROA');
        end
        [c,V] = tilqr(obj,obj.xG,obj.uG,Q,R);
        pp = feedback(obj.taylorApprox(0,obj.xG,obj.uG,3),c);
        options.method='levelSet';
        V=regionOfAttraction(pp,V,options);
      end
    end
    
    function [utraj,xtraj]=swingUpTrajectory(obj)
      x0 = zeros(14,1); 
      xf = double(obj.xG);
      tf0 = 4;
      
      N = 21;
      prog = DircolTrajectoryOptimization(obj,N,[2 6]);
      prog = prog.addStateConstraint(ConstantConstraint(x0),1);
      prog = prog.addStateConstraint(ConstantConstraint(xf),N);
      prog = prog.addRunningCost(@cost);
      prog = prog.addFinalCost(@finalCost);
      
      traj_init.x = PPTrajectory(foh([0,tf0],[double(x0),double(xf)]));
      
      for attempts=1:10
        tic
        [xtraj,utraj,z,F,info] = prog.solveTraj(tf0,traj_init);
        toc
        if info==1, break; end
      end

      function [g,dg] = cost(dt,x,u)
        R = 1;
        g = sum((R*u).*u,1);
        dg = [zeros(1,1+size(x,1)),2*u'*R];
        return;
        
        xd = repmat([pi;0;0;0],1,size(x,2));
        xerr = x-xd;
        xerr(1,:) = mod(xerr(1,:)+pi,2*pi)-pi;
        
        Q = diag([10,10,1,1]);
        R = 100;
        g = sum((Q*xerr).*xerr + (R*u).*u,1);
        
        if (nargout>1)
          dgddt = 0;
          dgdx = 2*xerr'*Q;
          dgdu = 2*u'*R;
          dg = [dgddt,dgdx,dgdu];
        end
      end
      
      function [h,dh] = finalCost(t,x)
        h = t;
        dh = [1,zeros(1,size(x,1))];
        return;
        
        xd = repmat([pi;0;0;0],1,size(x,2));
        xerr = x-xd;
        xerr(1,:) = mod(xerr(1,:)+pi,2*pi)-pi;
        
        Qf = 100*diag([10,10,1,1]);
        h = sum((Qf*xerr).*xerr,1);
        
        if (nargout>1)
          dh = [0, 2*xerr'*Qf];
        end
      end
    end
  end
  
end
