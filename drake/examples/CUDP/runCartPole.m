function runCartPole()

   % add all of the paths we need
   addpath('../CartPole');

   % get our plant
   p = CartPolePlant;
   v = CartPoleVisualizer(p);
      
   % set up traj_init and constants
   x0 = [0;0;0;0]; 
   xm = [1;pi/2;0;0];
   xf = [0;pi;0;0];
   tf0 = 4;
   N = 120;
   tourqueLimit = 30;
   Q = .1*eye(4);
   Qf = 1000*eye(4);
   R = 1*.01;
   SigmaScale = .01;
   tols = [[1e-2,1e-1];[1e-4,1e-2];[5e-7,5e-3]];
   traj_init.x = PPTrajectory(foh([0,tf0/2,tf0],[double(x0),double(xm),double(xf)]));
   
   % set up solvers
   options.debug = true;
   options.SigmaScale = SigmaScale;
   options.muMax = 1e30;
   options.phiFactor = 2;
   options.muFactor = 10;
   options.ILambdaDecay = 1;
   % penalty iLQR
   options.penaltyOnly = true;
   options.derivative_method = 1;
   options.tolConstr = tols(1,1);
   options.tolLambdaUpdate = tols(1,2);
   traj_opt1a = setUpTrajOpt(p,N,options,xf,tourqueLimit,Qf,Q,R);
   options.tolConstr = tols(2,1);
   options.tolLambdaUpdate = tols(2,2);
   traj_opt1b = setUpTrajOpt(p,N,options,xf,tourqueLimit,Qf,Q,R);
   options.tolConstr = tols(3,1);
   options.tolLambdaUpdate = tols(3,2);
   traj_opt1c = setUpTrajOpt(p,N,options,xf,tourqueLimit,Qf,Q,R);
   
   % penalty UDP
   options.tolConstr = tols(1,1);
   options.tolLambdaUpdate = tols(1,2);
   options.derivative_method = 2;
   traj_opt2a = setUpTrajOpt(p,N,options,xf,tourqueLimit,Qf,Q,R);
   options.tolConstr = tols(2,1);
   options.tolLambdaUpdate = tols(2,2);
   traj_opt2b = setUpTrajOpt(p,N,options,xf,tourqueLimit,Qf,Q,R);
   options.tolConstr = tols(3,1);
   options.tolLambdaUpdate = tols(3,2);
   traj_opt2c = setUpTrajOpt(p,N,options,xf,tourqueLimit,Qf,Q,R);
   
   % AL iLQR
   options.penaltyOnly = false;
   options.derivative_method = 1;
   options.tolConstr = tols(1,1);
   options.tolLambdaUpdate = tols(1,2);
   traj_opt3a = setUpTrajOpt(p,N,options,xf,tourqueLimit,Qf,Q,R);
   options.tolConstr = tols(2,1);
   options.tolLambdaUpdate = tols(2,2);
   traj_opt3b = setUpTrajOpt(p,N,options,xf,tourqueLimit,Qf,Q,R);
   options.tolConstr = tols(3,1);
   options.tolLambdaUpdate = tols(3,2);
   traj_opt3c = setUpTrajOpt(p,N,options,xf,tourqueLimit,Qf,Q,R);
      
   % AL UDP
   options.derivative_method = 2;
   options.tolConstr = tols(1,1);
   options.tolLambdaUpdate = tols(1,2);
   traj_opt4a = setUpTrajOpt(p,N,options,xf,tourqueLimit,Qf,Q,R);
   options.tolConstr = tols(2,1);
   options.tolLambdaUpdate = tols(2,2);
   traj_opt4b = setUpTrajOpt(p,N,options,xf,tourqueLimit,Qf,Q,R);
   options.tolConstr = tols(3,1);
   options.tolLambdaUpdate = tols(3,2);
   traj_opt4c = setUpTrajOpt(p,N,options,xf,tourqueLimit,Qf,Q,R);
      
   
   % solve
   % penalty iLQR
   tic;
   [x1a,u1a,d1a] = traj_opt1a.solveTraj(tf0,traj_init);
   t1a = toc;
   tic;
   [x1b,u1b,d1b] = traj_opt1b.solveTraj(tf0,traj_init);
   t1b = toc;
   tic;
   [x1c,u1c,d1c] = traj_opt1c.solveTraj(tf0,traj_init);
   t1c = toc;
   % penalty UDP
   tic;
   [x2a,u2a,d2a] = traj_opt2a.solveTraj(tf0,traj_init);
   t2a = toc;
   tic;
   [x2b,u2b,d2b] = traj_opt2b.solveTraj(tf0,traj_init);
   t2b = toc;
   tic;
   [x2c,u2c,d2c] = traj_opt2c.solveTraj(tf0,traj_init);
   t2c = toc;
   % AL iLQR
   tic;
   [x3a,u3a,d3a] = traj_opt3a.solveTraj(tf0,traj_init);
   t3a = toc;
   tic;
   [x3b,u3b,d3b] = traj_opt3b.solveTraj(tf0,traj_init);
   t3b = toc;
   tic;
   [x3c,u3c,d3c] = traj_opt3c.solveTraj(tf0,traj_init);
   t3c = toc;
   % AL UDP
   tic;
   [x4a,u4a,d4a] = traj_opt4a.solveTraj(tf0,traj_init);
   t4a = toc;
   tic;
   [x4b,u4b,d4b] = traj_opt4b.solveTraj(tf0,traj_init);
   t4b = toc;
   tic;
   [x4c,u4c,d4c] = traj_opt4c.solveTraj(tf0,traj_init);
   t4c = toc;
   
   keyboard;
   v.playback(x4c);
   
   function [g,dg,ddg] = cost(dt,x,u,xf,Q,R)
     g = (x-xf)'*Q*(x-xf) + u'*R*u;
     
     if (nargout>1)
       dg = [0, 2*(x-xf)'*Q,2*u'*R];
       ddg = blkdiag(0,2*Q,2*R);
     end
   end

   function [h,dh,ddh] = finalCost(tf,x,xf,Qf)
     h = (x-xf)'*Qf*(x-xf);
     if (nargout>1)
       dh = [0, 2*(x-xf)'*Qf];
       ddh = blkdiag(0,2*Qf);
     end
   end

   function [c,dc] = torqueLimit(k,u,tourqueLimit)
      c = [u - tourqueLimit; -u - tourqueLimit];
      dc = [1; -1];
   end

   function [c,dc] = xConstrEq(k,x,xf,N)
      if k < N
          c = [0 0 0 0]';
          dc = zeros(4,4);
      else
          c = [x-xf];
          dc = [eye(4)];
      end
   end

   function [traj_opt] = setUpTrajOpt(p,N,options,xf,tourqueLimit,Qf,Q,R)
       traj_opt = CUDPTrajectoryOptimization(p,N,[2 6],options);
       traj_opt = traj_opt.addInputConstraint(@(k,u)torqueLimit(k,u,tourqueLimit), 2);
       traj_opt = traj_opt.addStateConstraint(@(k,x)xConstrEq(k,x,xf,N), 4, 1);
       traj_opt = traj_opt.addRunningCost(@(t,x,u)cost(t,x,u,xf,Q,R));
       traj_opt = traj_opt.addFinalCost(@(t,x)finalCost(t,x,xf,Qf));
   end

end
