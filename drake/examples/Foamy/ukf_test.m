plant = FoamyPlant;

%trim conditions
[x0, u0] = findTrim(plant,7); %find trim conditions for level flight at 7 m/s
x0(3) = 1.5;

%Simulate at 50 Hz
tsamp = 0:.02:.4;
[tsamp, xsamp] = ode45(@(t,x) plant.dynamics(t,x,u0),tsamp,x0);
xsamp = xsamp';

%Generate sensor measurements
N = length(tsamp);
y = zeros(15,N);
for k = 1:N
    y(:,k) = foamy_sensors(xsamp(:,k),u0);
end

%Run UKF with offset on initial state
xhat = zeros(13,N);
xhat(:,1) = x0;
xhat(3,1) = xhat(3,1)-0.5;
Phat = zeros(12,12,N);
Phat(:,:,1) = blkdiag(1*eye(3), .01*eye(3), 0.1*eye(3), .01*eye(3));
Q = .001*eye(12);
R = .01*eye(15);
for k = 1:(N-1)
    [xhat(:,k+1), Phat(:,:,k+1)] = foamy_ukf(xhat(:,k),y(:,k+1),u0,Phat(:,:,k),Q,R,.02);
end


