%% Main file for the MPC control during the tracking of a target attitude.

clc
clear
close all

%% System parameters
param.mu = 398600.0 ;                   % Earth gravitational parameter
param.Re = 6378.0 ;                     % Earth radius
param.h = 500.0 ;                       % Orbit altitude
param.R = param.Re + param.h ;          % Orbit radius
param.n = sqrt(param.mu/(param.R^3)) ;  % Mean orbital motion

param.chaser_inertia = [230, 0, 0;
                        0, 240, 0;
                        0, 0, 25] ;
param.target_inertia = [0.5208, 0, 0;
                        0, 0.5208, 0;
                        0, 0, 0.6667] ;
param.wmax = 0.1 ; % [rad/s]

param.Jinv = inv(param.chaser_inertia);
                    
% Sampling time
Ts = 1 ; % [s]
param.Ts = Ts;

% Initial conditions:
o0          =   deg2rad([0; 0; 0]);
w_chaser    =   [randn(1); randn(1); randn(1)]*1e-3;

ic          =   [w_chaser;o0] ;

%% State space model + LQR
n = param.n ;

A = [0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0;
     1/4, 0, 0, 0, 0, 0;
     0, 1/4, 0, 0, 0, 0;
     0, 0, 1/4, 0, 0, 0];
     
J = param.chaser_inertia;
Bu = [1/J(1, 1), 0, 0;
     0, 1/J(2,2), 0;
     0, 0, 1/J(3, 3);
     0, 0, 0;
     0, 0, 0;
     0, 0, 0];
 
C = eye(6);

D = 0;
     
Q = .01*eye(6);
R = [1/(J(1, 1)^2), 0, 0;
     0, 1/(J(2, 2)^2), 0;
     0, 0, 1/(J(3, 3)^2)];

system          =   ss(A,Bu,C,D);

sys_d           =   c2d(system,Ts,'zoh');
[Ad,Bd,Cd,Dd]   =   ssdata(sys_d);

[K,S,e] = lqr(A,Bu,Q,R);

%% Open loop simulation
Nsim            =   100 ;

t               =   0:Ts:Nsim*Ts;
x0              =   ic ;
u0              =   [1e-4*sin(t/10); -1e-4*sin(t/10); 5e-5*sin(t/10)] ;

% Simulate with control input
[y,ti]          =   lsim(system,u0,t,x0);
figure(1),plot(ti,y(:,1:3),'Linewidth',1.5),grid on,hold on

[yd,ti]           =   lsim(sys_d,u0,t,x0);
plot(ti,yd(:,1:3),'--'),grid on,hold on

%% Open loop simulation for the nominal + output sensor noise system
Vd = 0.01*eye(6);   % disturbance covariance
Vn = .03*eye(size(C,1));   % noise covariance
param.Vd = Vd;
param.Vn = Vn;

Bf = [Bu, Vd, 0*Vn];

sysNoise = ss(A,Bf,C,[zeros(6,9),Vn]);

uDist = randn(6,size(t,2));
uNoise = randn(6,size(t,2));
u0 = [1e-4*sin(t/10); -1e-4*sin(t/10); 5e-5*sin(t/10)]  ;

uAug = [u0; Vd*Vd*uDist; Vn*Vn*uNoise];

[y_noise,ti,xn] = lsim(sysNoise, uAug, t,x0);
figure(2),plot(ti,y_noise(:,1:3)),grid on, hold on

figure(3),plot(ti,y_noise(:,4:6)),grid on, hold on

%% Open loop simulation for the nominal system + disturbance
sysDist = ss(A,Bf,eye(6),0);

[ydist,ti] = lsim(sysDist,uAug,t,x0);
figure(2),plot(ti,ydist(:,1:3),'k--'),grid on, hold on

figure(3),plot(ti,ydist(:,4:6),'k--'),grid on, hold on

%% Build Kalman Filter
[L,Pe,E] = lqe(A,Vd,C,Vd,Vn); % design kalman filter

sysKf = ss(A-L*C, [Bu,L], C, 0*[Bu,L]); % Kalman filter estimator

[yKF,ti] = lsim(sysKf,[u0;y_noise'],t,x0);
figure(2),plot(ti,yKF(:,1:3),'g','LineWidth',2)

xlabel('Time [s]'), ylabel('[rad/s]'), title('Open loop simulation with KF')
legend('Noisy Measurement', '','','Model Measurment', '','','Filtered Measurement')

figure(3),plot(ti,yKF(:,4:6),'g','LineWidth',2)

xlabel('Time [s]'), ylabel('[rad]'), title('Open loop simulation with KF')
legend('Noisy Measurement', '','','Model Measurment', '','','Filtered Measurement')

%% Reference attitude (qauternion + angular velocity)
N_MPC       =   50;
N           =   Nsim+N_MPC ;
Tt          =   0:Ts:Ts*(N-1);

qt          =   [0; 0; 0];
w_target    =   [0.0; 0.0; 0.01];
target_ic   =   [w_target;qt];

x_ref        = zeros(size(ic,1),N);
x_ref(1:6,1) = target_ic;
for i =2:N
    x_ref(1:6,i) = targetode(x_ref(1:6, i-1),[0;0;0],param);
end

% Plot target mrp evolution
figure(3),subplot(4,1,1),plot(Tt,x_ref(4,:)), grid on, hold on, title('o1 ref')
figure(3),subplot(4,1,2),plot(Tt,x_ref(5,:)), grid on, hold on, title('o2 ref')
figure(3),subplot(4,1,3),plot(Tt,x_ref(6,:)), grid on, hold on, title('o3 ref')

% Plot target angular velocity evolution
figure(4),subplot(3,1,1),plot(Tt,x_ref(1,:)), grid on, hold on, title('wx ref')
figure(4),subplot(3,1,2),plot(Tt,x_ref(2,:)), grid on, hold on, title('wy ref')
figure(4),subplot(3,1,3),plot(Tt,x_ref(3,:)), grid on, hold on, title('wz ref')

close all

%% LQR Control Problem

xsim            =   zeros(size(ic,1),Nsim);
xsim(:,1)       =   x0;
xdsim           =   zeros(size(ic,1),Nsim);
xdsim(:,1)      =   x0;
ydsim           =   zeros(size(ic,1),Nsim);
ydsim(:,1)      =   C*x0;
xHatsim         =   zeros(size(ic,1),Nsim);
xHatsim(:,1)    =   x0;
tsim            =   0:Ts:(Nsim-1)*Ts;

ulqr        =   zeros(Nsim,3);
uMax        =   0.01;
uMin        =   -0.01;

d           =   [0;0;0] ;

for ind_sim=2:Nsim
    u                   =   - K * (xHatsim(:,ind_sim-1) - x_ref(:,ind_sim-1));
    u(u>uMax)           =   uMax;
    u(u<uMin)           =   uMin;

%     d                   =   [normrnd(0, 1e-3); normrnd(0, 1e-3); normrnd(0, 1e-3)] ;

    % Integrate nominal dynamics
    xout    = linearAdcsModel(xsim(:, ind_sim-1),u,d,param) ;
    xsim(:,ind_sim) = xout ;

    % Integrate disturbed dynamic measurements
    uAug = [u; Vd*Vd*randn(6,1); Vn*Vn*randn(6,1)];  % !!!!

    [tout1, xdout]   = ode45(@(t,x)linearAdcsModelDist(t,x,uAug,sysNoise), [0 Ts], xsim(:, ind_sim-1)) ;
    xdsim(:,ind_sim) = sysNoise.C*xdout(end,:)';
    ydsim(:,ind_sim) = sysNoise.C * xdsim(:,ind_sim);

    % Estimate the state
    [tout3, xkout]      = ode45(@(t,x)linearAdcsModelKalman(t,x,u,ydsim(:,ind_sim),sysKf), [0 Ts], xHatsim(:, ind_sim-1)) ;
    xHatsim(:,ind_sim)  = xkout(end,:)';

    ulqr(ind_sim-1,:)   =   u;
end

% Angular velocity vs reference
figure(5),subplot(3,1,1),plot(tsim,ydsim(1,:)),grid on, hold on, title('wx [rad/s]')
figure(5),subplot(3,1,2),plot(tsim,ydsim(2,:)),grid on, hold on, title('wy [rad/s]')
figure(5),subplot(3,1,3),plot(tsim,ydsim(3,:)),grid on, hold on, title('wz [rad/s]')

figure(5),subplot(3,1,1),plot(tsim,xHatsim(1,:)),grid on, hold on, title('wx [rad/s]')
figure(5),subplot(3,1,2),plot(tsim,xHatsim(2,:)),grid on, hold on, title('wy [rad/s]')
figure(5),subplot(3,1,3),plot(tsim,xHatsim(3,:)),grid on, hold on, title('wz [rad/s]')

figure(5),subplot(3,1,1),plot(tsim,x_ref(1,1:end-N_MPC),'r--'),grid on, hold on
figure(5),subplot(3,1,2),plot(tsim,x_ref(2,1:end-N_MPC),'r--'),grid on, hold on
figure(5),subplot(3,1,3),plot(tsim,x_ref(3,1:end-N_MPC),'r--'),grid on, hold on

legend('Sensor','KF','Ref')

% Angle vs reference
figure(6),subplot(3,1,1),plot(tsim,ydsim(4,:)),grid on, hold on, title('x [deg]')
figure(6),subplot(3,1,2),plot(tsim,ydsim(5,:)),grid on, hold on, title('y [deg]')
figure(6),subplot(3,1,3),plot(tsim,ydsim(6,:)),grid on, hold on, title('z [deg]')

figure(6),subplot(3,1,1),plot(tsim,xHatsim(4,:)),grid on, hold on, title('x [deg]')
figure(6),subplot(3,1,2),plot(tsim,xHatsim(5,:)),grid on, hold on, title('y [deg]')
figure(6),subplot(3,1,3),plot(tsim,xHatsim(6,:)),grid on, hold on, title('z [deg]')

figure(6),subplot(3,1,1),plot(tsim,x_ref(4,1:end-N_MPC),'r--'),grid on, hold on
figure(6),subplot(3,1,2),plot(tsim,x_ref(5,1:end-N_MPC),'r--'),grid on, hold on
figure(6),subplot(3,1,3),plot(tsim,x_ref(6,1:end-N_MPC),'r--'),grid on, hold on

% close all

figure(7),subplot(2,1,1),plot(tsim,ulqr(:,1)),grid on, hold on, title('Input acceleration [km/s^2] with LQR vs MPC')
figure(7),subplot(2,1,1),plot(tsim,ulqr(:,2)),grid on, hold on
figure(7),subplot(2,1,1),plot(tsim,ulqr(:,3)),grid on, hold on
