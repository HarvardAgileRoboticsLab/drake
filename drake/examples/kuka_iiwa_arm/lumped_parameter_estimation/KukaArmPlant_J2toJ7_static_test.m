classdef KukaArmPlant_J2toJ7_static_test < Manipulator 
  
  properties
    % parameters from KUKA CAD Model with emperical tunning
        
    l1x = 0, l1y = 0, l1z = 0.1575;
    l2x = 0, l2y = 0, l2z = 0.2025;
    l3x = 0, l3y = 0.2045, l3z = 0;
    l4x = 0, l4y = 0, l4z = 0.2155;
    l5x = 0, l5y = 0.1845, l5z = 0;
    l6x = 0, l6y = 0, l6z = 0.2155;
    l7x = 0,l7y = 0.081, l7z = 0;
    
    m1 = 5.76; 
    m2 = 6.35; 
    m3 = 3.5;
    m4 = 3.5; 
    m5 = 3.5; 
    m6 = 1.8; m7 = 1.18;
    g = 9.81;
    
    % viscous coefficients of positive and negative directions
    bv1_positive=0.1392;  bv2_positive=8.3008; bv3_positive=0.0108;  bv4_positive=1.4523; bv5_positive=0.0178;  bv6_positive=0.1199; bv7_positive= -0.0573;
    bv1_negative=0.2439;  bv2_negative=8.6682; bv3_negative=-0.0106;  bv4_negative=1.6164; bv5_negative=0.0104;  bv6_negative=0.0961; bv7_negative= -0.0319;
    % coulomb coefficients of positive and negative directions
    bc1_positive=1.7809;  bc2_positive= 1.9103; bc3_positive= -0.0445;  bc4_positive=0.2538; bc5_positive=0.1151;  bc6_positive=0.0534; bc7_positive=0.4934;
    bc1_negative= -0.3484;  bc2_negative= 0.1029; bc3_negative= -0.3254;  bc4_negative=-0.4345; bc5_negative=0.0809;  bc6_negative=-0.0881; bc7_negative=-0.1781;    
    
    c1x = 0, c1y = -0.03, c1z = 0.12;
    c2x = 0.0003, c2y = 0.059, c2z = 0.042;
    c3x = 0, c3y = 0.03, c3z = 0.13;
    c4x = 0, c4y = 0.067, c4z = 0.034;
    c5x = 0.0001, c5y = 0.021, c5z = 0.076;
    c6x = 0, c6z = 0.0004;
    c7x = 0, c7y = 0;

    % update with identified parameters
    c7z = 0.0176;%0.02;
    c6y = 0.00055;%0.0006 
    
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
    function obj = KukaArmPlant_J2toJ7_static_test
      obj = obj@Manipulator(6,6);
      obj = setInputLimits(obj,-10,10);

      obj = setInputFrame(obj,CoordinateFrame('KukaArmInput',6,'u',{'tau2','tau3','tau4','tau5','tau6','tau7'}));
      obj = setStateFrame(obj,CoordinateFrame('KukaArmState',12,'x',{'theta2','theta3','theta4','theta5','theta6','theta7', ...
                                                'theta2dot','theta3dot','theta4dot','theta5dot','theta6dot','theta7dot'}));
      obj = setOutputFrame(obj,obj.getStateFrame);
      
      obj.xG = Point(obj.getStateFrame,[pi/2;pi/2;pi/2;pi/2;pi/2;pi/2;0;0;0;0;0;0]);%[To Be Checked]
      obj.uG = Point(obj.getInputFrame,[0;0;0;0;0;0]);

      obj = setParamFrame(obj,CoordinateFrame('KukaArmParams',10,'p',...
        {'m2','m3','m4','m5','c2x', 'c2y','c3y','c3z','c4y','c4z'})); %'m2','m3','m4','m1','m2','m3',,'m6','m7' 'c2x', 'c2y', 'c2z' 
      obj = setParamLimits(obj,zeros(obj.getParamFrame.dim,1));
    end

    function [C,B] = manipulatorDynamics(obj,q,qd)
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
      
      tic
      c = cos(q(1:6,:));  s = sin(q(1:6,:));%[To be Checked]
      toc

      % Inertia matrix and Coriolis and centrifugal vector are ignored
      % since this is only static test, qddot = 0, qdot = 0.
      tic
      %H = getKukaArmInertiaMatrix_J2toJ7(obj,q,c,s);
      toc
      
      tic
      %C = getKukaArmCoriolisCentrifugalVector_J2toJ7(obj,q,qd,c,s);
      toc
      tic
      G = getKukaArmGravityVector_J2toJ7(obj,q,c,s);
      toc

      b(1,1) = bv2_positive * qd(1);
      b(2,1) = bv3_positive * qd(2);
      b(3,1) = bv4_positive * qd(3);
      b(4,1) = bv5_positive * qd(4);
      b(5,1) = bv6_positive * qd(5);
      b(6,1) = bv7_positive * qd(6);
      
      % accumate total C and add a damping term:
      C = G;
      B = eye(6);
    end
    
    % to be updated
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
      
      q = x(1:7); qd = x(8:14);
      c = cos(q(1:7,:));  s = sin(q(1:7,:));
            
      T = getKineticEnergy(obj,q,qd);
      U = getPotentialEnergy(obj,q);      
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
      x = .1*randn(10,1);
    end
    
    function n = getNumPositions(obj)
      n = 5;
    end
    
    function n = getNumVelocities(obj)
      n = 5;
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