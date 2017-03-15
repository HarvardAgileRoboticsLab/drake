classdef KukaArmPlant_J5toJ7_dynamic_test < Manipulator 
  
  properties
    % parameters from KUKA CAD Model with emperical tunning
    
    l5x = 0, l5y = 0.1845, l5z = 0;
    l6x = 0, l6y = 0, l6z = 0.2155;
    l7x = 0,l7y = 0.081, l7z = 0;
    
    m5 = 3.5; m6 = 1.8; m7 = 1.2;
    g = 9.81;
    
    % viscous coefficients of positive and negative directions
    bv1_positive=0.1392;  bv2_positive=8.3008; bv3_positive=0.0108;  bv4_positive=1.4523; bv5_positive=0.0178;  bv6_positive=0.1199; bv7_positive= -0.0573;
    bv1_negative=0.2439;  bv2_negative=8.6682; bv3_negative=-0.0106;  bv4_negative=1.6164; bv5_negative=0.0104;  bv6_negative=0.0961; bv7_negative= -0.0319;
    % coulomb coefficients of positive and negative directions
    bc1_positive=1.7809;  bc2_positive= 1.9103; bc3_positive= -0.0445;  bc4_positive=0.2538; bc5_positive=0.1151;  bc6_positive=0.0534; bc7_positive=0.4934;
    bc1_negative= -0.3484;  bc2_negative= 0.1029; bc3_negative= -0.3254;  bc4_negative=-0.4345; bc5_negative=0.0809;  bc6_negative=-0.0881; bc7_negative=-0.1781;

    c5x = 0.0001, c5y = 0.021, c5z = 0.076;
    c6x = 0, c6y = 0.0006, c6z = 0.0004;
    c7x = 0, c7y = 0, c7z = 0.02;
    
    I5xx= 0.02, I5xy= 0, I5xz= 0, I5yy= 0.018, I5yz= 0.00, I5zz= 0.005;
    I6xx= 0.005, I6xy= 0, I6xz= 0, I6yy= 0.0036, I6yz= 0.000, I6zz= 0.0047;
    I7xx= 0.001, I7xy= 0, I7xz= 0, I7yy= 0.001, I7yz= 0, I7zz= 0.001;
    % by a rough measurement (treat the 7th joint as a cylinder, I7zz is roughtly a number between 0.0004 and 0.0014)

    xG
    uG
  end
  
  methods
    function obj = KukaArmPlant_J5toJ7_dynamic_test
      obj = obj@Manipulator(3,3);
      obj = setInputLimits(obj,-10,10);

      obj = setInputFrame(obj,CoordinateFrame('KukaArmInput',3,'u',{'tau5','tau6','tau7'}));
      obj = setStateFrame(obj,CoordinateFrame('KukaArmState',6,'x',{'theta5','theta6','theta7', ...
                                                'theta5dot','theta6dot','theta7dot'}));
      obj = setOutputFrame(obj,obj.getStateFrame);
      
      obj.xG = Point(obj.getStateFrame,[0;0;0;0;0;0]);%[To Be Checked]
      obj.uG = Point(obj.getInputFrame,[0;0;0]);

      obj = setParamFrame(obj,CoordinateFrame('KukaArmParams',7,'p',...
        {'I5zz','I6xx','I6yy','I6zz','I7xx','I7yy','I7zz'}));
      obj = setParamLimits(obj,zeros(obj.getParamFrame.dim,1));
    end

    function [H,C,B] = manipulatorDynamics(obj,q,qd)
      % keep it readable:
      m5=obj.m5; m6=obj.m6; m7=obj.m7;
      l5x=obj.l5x; l5y=obj.l5y; l5z=obj.l5z;
      l6x=obj.l6x; l6y=obj.l6y; l6z=obj.l6z;
      l7x=obj.l7x; l7y=obj.l7y; l7z=obj.l7z;
      
      % link inertia (with repsect to joint local coordinate)
      I5xx= obj.I5xx; I5xy= obj.I5xy; I5xz= obj.I5xz; I5yy= obj.I5yy; I5yz= obj.I5yz; I5zz= obj.I5zz;
      I6xx= obj.I6xx; I6xy= obj.I6xy; I6xz= obj.I6xz; I6yy= obj.I6yy; I6yz= obj.I6yz; I6zz= obj.I6zz;
      I7xx= obj.I7xx; I7xy= obj.I7xy; I7xz= obj.I7xz; I7yy= obj.I7yy; I7yz= obj.I7yz; I7zz= obj.I7zz;

      g=obj.g; 
      c5x = obj.c5x; c5y = obj.c5y; c5z = obj.c5z;
      c6x = obj.c6x; c6y = obj.c6y; c6z = obj.c6z;
      c7x = obj.c7x; c7y = obj.c7y; c7z = obj.c7z;
      
      bv1_positive=obj.bv1_positive; bv2_positive=obj.bv2_positive; bv3_positive=obj.bv3_positive; bv4_positive=obj.bv4_positive; bv5_positive=obj.bv5_positive; bv6_positive=obj.bv6_positive; bv7_positive=obj.bv7_positive;
      bv1_negative=obj.bv1_negative; bv2_negative=obj.bv2_negative; bv3_negative=obj.bv3_negative; bv4_negative=obj.bv4_negative; bv5_negative=obj.bv5_negative; bv6_negative=obj.bv6_negative; bv7_negative=obj.bv7_negative;
      bc1_positive=obj.bc1_positive; bc2_positive=obj.bc2_positive; bc3_positive=obj.bc3_positive; bc4_positive=obj.bc4_positive; bc5_positive=obj.bc5_positive; bc6_positive=obj.bc6_positive; bc7_positive=obj.bc7_positive;
      bc1_negative=obj.bc1_negative; bc2_negative=obj.bc2_negative; bc3_negative=obj.bc3_negative; bc4_negative=obj.bc4_negative; bc5_negative=obj.bc5_negative; bc6_negative=obj.bc6_negative; bc7_negative=obj.bc7_negative;
            
      q5 = q(1);
      q6 = q(2);
      q7 = q(3);
      
      % joint angular velocity data
      dq5 = qd(1);
      dq6 = qd(2);
      dq7 = qd(3);      
      
      tic      
      h11 = I5zz + m6*(c6y^2*sin(q6)^2 + c6z^2) + m5*(c5x^2 + c5y^2) + I6yy*cos(q6)^2 + I7zz*cos(q6)^2 + I6xx*sin(q6)^2 + I7xx*cos(q7)^2*sin(q6)^2 + I7yy*sin(q6)^2*sin(q7)^2 + (m7*sin(q6)^2*(1000*c7z + 81)^2)/1000000;
      h21 = (2268949521066275*cos(q6))/9223372036854775808 - I7xx*cos(q7)*sin(q6)*sin(q7) + I7yy*cos(q7)*sin(q6)*sin(q7) - c6y*c6z*m6*cos(q6);
      h22 = m6*c6y^2 + m7*c7z^2 + (81*m7*c7z)/500 + I7xx/2 + I7yy/2 + I6zz + (6561*m7)/1000000 - (I7xx*cos(2*q7))/2 + (I7yy*cos(2*q7))/2;
      h31 = I7zz*cos(q6);
      h32 = 0;
      h33 = I7zz;
      h12 = h21; h13 = h31; h23 = h32;
      H = [h11, h12, h13;h21, h22, h23;h31, h32, h33];
      toc
      
      tic
      C1 = (6561*dq5*dq6*m7*sin(2*q6))/1000000 - I7zz*dq6*dq7*sin(q6) - (2268949521066275*dq5*dq6*sin(q5))/9223372036854775808 + c6y*c6z*dq6^2*m6*sin(q6) + (81*c7z*dq5*dq6*m7*sin(2*q6))/500 + c6y^2*dq5*dq6*m6*sin(2*q6) + c7z^2*dq5*dq6*m7*sin(2*q6);
      C2 = dq5*((I7xx*dq6*sin(2*q5))/2 - (I7yy*dq6*sin(2*q5))/2 + I7yy*dq7*cos(q5)^2*sin(q6) + I7xx*dq7*sin(q5)^2*sin(q6)) - (dq5^2*m7*sin(2*q6)*(1000*c7z + 81)^2)/2000000 - (c6y^2*dq5^2*m6*sin(2*q6))/2 + (dq5*dq6*sin(2*q5)*(I6xx - I6yy))/2 - dq6*dq7*cos(q5)*cos(q6)*sin(q5)*(I7xx - I7yy);
      C3 = dq6*dq7*cos(q6)*sin(q6)*(I7xx - I7zz - I7xx*sin(q5)^2 + I7yy*sin(q5)^2) - dq5*sin(q6)*(I7xx*dq6*cos(q5)^2 + I7yy*dq6*sin(q5)^2 + I7xx*dq7*cos(q5)*sin(q5)*sin(q6) - I7yy*dq7*cos(q5)*sin(q5)*sin(q6));
      C = [C1;C2;C3];
      toc
      tic
      G1 = 0;
      G2 = -(g*sin(q6)*(81*m7 + 1000*c6y*m6 + 1000*c7z*m7))/1000; % -(3*g*m6*sin(q6))/5000 - (101*g*m7*sin(q6))/1000;
      G3 = 0;
      G = [G1;G2;G3];
      toc
      % Coulomb and viscous friction
      
      b1 = bv5_positive * dq5;% + bc5_positive;
      b2 = bv6_positive * dq6;% + bc6_positive;      
      b3 = bv7_positive * dq7;% + bc7_positive;
      b = [b1;b2;b3];
      C = C + G;
      B = eye(3);
    end
    
    % to be updated
    function [T,U] = energy(obj,x)
      m5=obj.m5; m6=obj.m6; m7=obj.m7; 
      
      l5x=obj.l5x; l5y=obj.l5y; l5z=obj.l5z;
      l6x=obj.l6x; l6y=obj.l6y; l6z=obj.l6z;
      l7x=obj.l7x; l7y=obj.l7y; l7z=obj.l7z;
      
      % link inertia (with repsect to joint local coordinate)
      I5xx= obj.I5xx; I5xy= obj.I5xy; I5xz= obj.I5xz; I5yy= obj.I5yy; I5yz= obj.I5yz; I5zz= obj.I5zz;
      I6xx= obj.I6xx; I6xy= obj.I6xy; I6xz= obj.I6xz; I6yy= obj.I6yy; I6yz= obj.I6yz; I6zz= obj.I6zz;
      I7xx= obj.I7xx; I7xy= obj.I7xy; I7xz= obj.I7xz; I7yy= obj.I7yy; I7yz= obj.I7yz; I7zz= obj.I7zz;

      g=obj.g; 
      c5x = obj.c5x; c5y = obj.c5y; c5z = obj.c5z;
      c6x = obj.c6x; c6y = obj.c6y; c6z = obj.c6z;
      c7x = obj.c7x; c7y = obj.c7y; c7z = obj.c7z;
      
      q = x(1:3); qd = x(4:6);
      
      q5 = q(1);
      q6 = q(2);
      q7 = q(3);

      % joint angular velocity data
      dq5 = qd(1);
      dq6 = qd(2);
      dq7 = qd(3);
      
      T = (I6xx*dq5^2)/2 + (I7xx*dq5^2)/4 + (I7xx*dq6^2)/4 + (I7yy*dq5^2)/4 + (I7yy*dq6^2)/4 + (I5zz*dq5^2)/2 + (I6zz*dq6^2)/2 + (I7zz*dq7^2)/2 + (44101*dq5^2*m5)/200000000 + (13*dq5^2*m6)/50000000 + (10201*dq5^2*m7)/2000000 + (9*dq6^2*m6)/50000000 + (10201*dq6^2*m7)/2000000 + (I7xx*dq5^2*cos(2*q7))/4 - (I7xx*dq6^2*cos(2*q7))/4 - (I7yy*dq5^2*cos(2*q7))/4 + (I7yy*dq6^2*cos(2*q7))/4 - (I6xx*dq5^2*cos(q6)^2)/2 - (I7xx*dq5^2*cos(q6)^2)/4 + (I6yy*dq5^2*cos(q6)^2)/2 - (I7yy*dq5^2*cos(q6)^2)/4 + (I7zz*dq5^2*cos(q6)^2)/2 - (9*dq5^2*m6*cos(q6)^2)/50000000 - (10201*dq5^2*m7*cos(q6)^2)/2000000 + (2268949521066275*dq5*dq6*cos(q6))/9223372036854775808 + I7zz*dq5*dq7*cos(q6) - (3*dq5*dq6*m6*cos(q6))/12500000 - (I7xx*dq5^2*cos(2*q7)*cos(q6)^2)/4 + (I7yy*dq5^2*cos(2*q7)*cos(q6)^2)/4 - (I7xx*dq5*dq6*sin(2*q7)*sin(q6))/2 + (I7yy*dq5*dq6*sin(2*q7)*sin(q6))/2;
      U = (19*g*m5)/250 + g*m7*((101*cos(q6))/1000 + 431/2000) + g*m6*((3*cos(q6))/5000 + 431/2000);
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
      x = .1*randn(6,1);
    end
    
    function n = getNumPositions(obj)
      n = 3;
    end
    
    function n = getNumVelocities(obj)
      n = 3;
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
