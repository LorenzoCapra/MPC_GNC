%% Main file for the MPC control during the tracking of a target FULL pose.

clc
clear
close all

%% System parameters
param.mu = 398600.0 ;                   % Earth gravitational parameter
param.Re = 6378.0 ;                     % Earth radius
param.h = 500.0 ;                       % Orbit altitude
param.R = param.Re + param.h ;          % Orbit radius
param.n = sqrt(param.mu/(param.R^3)) ;  % Mean orbital motion
param.a = 7000 ;                        % Target semi-major axis
param.chaser_inertia = [0.872, 0, 0;
                        0, 0.115, 0;
                        0, 0, 0.797] ;  % Chaser inertia matrix
param.target_inertia = [0.5208, 0, 0;
                        0, 0.5208, 0;
                        0, 0, 0.6667] ; % Target inertia matrix
                   
param.wmax = 0.1 ;                      % Max angualar velocity [rad/s]

% Sampling time
Ts = 1 ; % [s]

% Initial conditions:
alfa        =   [0; 0; 0];
w_target    =   [0.0; 0.0; 0.05];
target_ic   =   [alfa; w_target];

r_mag0      =   0.003;
r0          =   [r_mag0; 0; 0];
v0          =   cross(w_target, r0);

q0          =   [0; 0; 0; 1];
w_chaser    =   [0; 0; 0];

ic          =   [r0; v0; q0; w_chaser] ;

%% Open loop simulation
Nsim            =   30 ;
x0              =   ic ;
u               =   [normrnd(0, 1e-3); normrnd(0, 1e-3); normrnd(0, 1e-3);
                     normrnd(0, 1e-3); normrnd(0, 1e-3); normrnd(0, 1e-3)] ;
d               =   [normrnd(0, 1e-4); normrnd(0, 1e-4); normrnd(0, 1e-4);
                     normrnd(0, 1e-4); normrnd(0, 1e-4); normrnd(0, 1e-4)];

xsim            =   zeros(size(ic,1), Nsim) ;
xsim(:,1)       =   x0 ;
tsim            =   0:Ts:(Nsim-1)*Ts ;

for ind_sim=2:Nsim
    [tout, xout]     = ode45(@(t,x)fulltrackmodel(t,x,u,d,param), [0 Ts], xsim(:,ind_sim-1)) ;
    xsim(:,ind_sim)  = xout(end,:)' ;
end

figure(1),subplot(3,1,1),plot(tsim,xsim(1,:)),grid on, hold on, title('x [km]')
figure(1),subplot(3,1,2),plot(tsim,xsim(2,:)),grid on, hold on, title('y [km]')
figure(1),subplot(3,1,3),plot(tsim,xsim(3,:)),grid on, hold on, title('z [km]')

figure(2),subplot(4,1,1),plot(tsim,xsim(7,:)),grid on, hold on, title('q1')
figure(2),subplot(4,1,2),plot(tsim,xsim(8,:)),grid on, hold on, title('q2')
figure(2),subplot(4,1,3),plot(tsim,xsim(9,:)),grid on, hold on, title('q3')
figure(2),subplot(4,1,4),plot(tsim,xsim(10,:)),grid on, hold on, title('q4')

figure(3),subplot(3,1,1),plot(tsim,xsim(11,:)),grid on, hold on, title('wc_x [rad/s]')
figure(3),subplot(3,1,2),plot(tsim,xsim(12,:)),grid on, hold on, title('wc_y [rad/s]')
figure(3),subplot(3,1,3),plot(tsim,xsim(13,:)),grid on, hold on, title('wc_z [rad/s]')

figure(4),plot3(xsim(1,:),xsim(2,:),xsim(3,:),'LineWidth',1.5),grid on,hold on,
title('Relative Motion Chaser-Target')
xlabel('x [km]'), ylabel('y [km]'), zlabel('z [km]')
figure(4),scatter3(0,0,0, 'filled'),grid on

%% Finite Horizon Optimal Control Problem
x0          =   ic ;
N           =   40 ;
Q           =   1 ;
R           =   1e-5 ;
umax        =   [1e-2 1e-2 1e-2 1e-2 1e-2 1e-2] ;
d           =   [0;0;0;0;0;0];

% ---------- Reference trajectory (or point) --------------------------------------------------------- %
N_MPC       =   10;
x_r         =   zeros(12,N+N_MPC);

% Reference position/velocity
phi         =   linspace(0,w_target(3)*N,N+N_MPC);
%r_mag      =   linspace(0.01, 0.003, N);
r_ref       =   zeros(3,N+N_MPC);
v_ref       =   zeros(3,N+N_MPC);
for i = 1:N+N_MPC
    r_ref(:,i)      = [r_mag0 * cos(phi(i)),r_mag0 * sin(phi(i)), 0];
    v_ref(:,i)      = cross(w_target,r_ref(:,i));
    
    x_r(1:6,i)      = [r_ref(:,i); v_ref(:,i)];
end

% Reference attitude/angular velocity
Tt              = 0:Ts:Ts*(N-1+N_MPC);
[~,xtout]       = ode45(@(t,x)targetode(t,x,param),Tt,target_ic);
att_ref         = xtout';

x_r(7:end,:)  = att_ref(1:6,:);

figure(5),subplot(4,1,1),plot3(r_ref(1,:),r_ref(2,:),r_ref(3,:)),grid on,hold on,title('Reference trajectory')
figure(5),subplot(4,1,1),scatter3(0,0,0,'filled')
figure(5),subplot(4,1,2),plot3(v_ref(1,:),v_ref(2,:),v_ref(3,:)),grid on,title('Reference velocity')
figure(5),subplot(4,2,1),plot(Tt,att_ref(1,:)),grid on,title('Reference Angles')
figure(5),subplot(4,2,1),plot(Tt,att_ref(2,:)),grid on,title('Reference Angles')
figure(5),subplot(4,2,1),plot(Tt,att_ref(3,:)),grid on,title('Reference Angles')
figure(5),subplot(4,2,2),plot(Tt,att_ref(4,:)),grid on,title('Reference Angular Velocity')
figure(5),subplot(4,2,2),plot(Tt,att_ref(5,:)),grid on,title('Reference Angular Velocity')
figure(5),subplot(4,2,2),plot(Tt,att_ref(6,:)),grid on,title('Reference Angular Velocity')

% Transform euler angles to quaternion
x_ref         =   zeros(13,N+N_MPC);
for i = 1:N+N_MPC
    x_ref_quat  =   DCM_quat(Angles321_DCM(att_ref(1:3,i))) ;
    x_ref(:,i)  =   [x_r(1:6,i); x_ref_quat; x_r(10:12,i)];
end

% ---------------------------------------------------------------------------------------------------- %

%% Set the options and solve the optimization problem in open loop
options     =   optimset('Display','Iter','MaxFunEvals',1e4,'Algorithm','active-set');  

ustar   = fmincon(@(u)fulltrack_cost_fun(x0,u,d,N,Q,R,x_ref,1,Ts,param),zeros(N,6),...
                [],[],[],[],[],[],@(u)fulltrack_constr_fun(x0,u,d,N,umax,Ts,param),options);

%% Open loop simulation - FHOCP solution
Nsim            =   N;
d               =   [0;0;0;0;0;0];

xsim            =   zeros(size(ic,1),Nsim);
xsim(:,1)       =   x0;
tsim            =   0:Ts:(Nsim-1)*Ts;

for ind_sim=2:Nsim
    u               = ustar(ind_sim-1,:)';
    [tout,xout]     = ode45(@(t,x)fulltrackmodel(t,x,u,d,param),[0 Ts],xsim(:,ind_sim-1));
    xsim(:,ind_sim) = xout(end,:)';
end

% Position vs reference
figure(6),subplot(3,1,1),plot(tsim,xsim(1,:)),grid on, hold on, title('x [km]')
figure(6),subplot(3,1,2),plot(tsim,xsim(2,:)),grid on, hold on, title('y [km]')
figure(6),subplot(3,1,3),plot(tsim,xsim(3,:)),grid on, hold on, title('z [km]')

figure(6),subplot(3,1,1),plot(tsim,x_ref(1,1:end-N_MPC),'r--'),grid on, hold on
figure(6),subplot(3,1,2),plot(tsim,x_ref(2,1:end-N_MPC),'r--'),grid on, hold on
figure(6),subplot(3,1,3),plot(tsim,x_ref(3,1:end-N_MPC),'r--'),grid on, hold on

% Velocity vs reference
figure(7),subplot(3,1,1),plot(tsim,xsim(4,:)),grid on, hold on, title('vx [km/s]')
figure(7),subplot(3,1,2),plot(tsim,xsim(5,:)),grid on, hold on, title('vy [km/s]')
figure(7),subplot(3,1,3),plot(tsim,xsim(6,:)),grid on, hold on, title('vz [km/s]')

figure(7),subplot(3,1,1),plot(tsim,x_ref(4,1:end-N_MPC),'r--'),grid on, hold on
figure(7),subplot(3,1,2),plot(tsim,x_ref(5,1:end-N_MPC),'r--'),grid on, hold on
figure(7),subplot(3,1,3),plot(tsim,x_ref(6,1:end-N_MPC),'r--'),grid on, hold on

% Angles vs reference
figure(8),subplot(3,1,1),plot(tsim,xsim(7,:)),grid on, hold on, title('alfa-x [rad]')
figure(8),subplot(3,1,2),plot(tsim,xsim(8,:)),grid on, hold on, title('alfa-y [rad]')
figure(8),subplot(3,1,3),plot(tsim,xsim(9,:)),grid on, hold on, title('alfa-z [rad]')

figure(8),subplot(3,1,1),plot(tsim,x_ref(7,1:end-N_MPC),'r--'),grid on, hold on
figure(8),subplot(3,1,2),plot(tsim,x_ref(8,1:end-N_MPC),'r--'),grid on, hold on
figure(8),subplot(3,1,3),plot(tsim,x_ref(9,1:end-N_MPC),'r--'),grid on, hold on

% Angular Velocity vs reference
figure(9),subplot(3,1,1),plot(tsim,xsim(10,:)),grid on, hold on, title('wx [rad/s]')
figure(9),subplot(3,1,2),plot(tsim,xsim(11,:)),grid on, hold on, title('wy [rad/s]')
figure(9),subplot(3,1,3),plot(tsim,xsim(12,:)),grid on, hold on, title('wz [rad/s]')

figure(9),subplot(3,1,1),plot(tsim,x_ref(10,1:end-N_MPC),'r--'),grid on, hold on
figure(9),subplot(3,1,2),plot(tsim,x_ref(11,1:end-N_MPC),'r--'),grid on, hold on
figure(9),subplot(3,1,3),plot(tsim,x_ref(12,1:end-N_MPC),'r--'),grid on, hold on

% Plot control effort
figure(666),subplot(3,1,1),plot(tsim,ustar(:,1)),grid on, hold on,
figure(666),subplot(3,1,2),plot(tsim,ustar(:,2)),grid on, hold on,
figure(666),subplot(3,1,3),plot(tsim,ustar(:,3)),grid on, hold on

%% MPC (Close loop simulation w. receding horizon)
Nsim            =   N;
x0              =   ic ;
d               =   [0;0;0;0;0;0];

xMPC            =   zeros(size(ic,1),Nsim);
uMPC            =   zeros(Nsim,6);
xMPC(:,1)       =   x0;
tsim            =   0:Ts:(Nsim-1)*Ts;
ustar           =   [zeros(N,6);zeros(1,6)];

for ind_sim=2:Nsim
    ustar               =   fmincon(@(u)fulltrack_cost_fun(xMPC(:,ind_sim-1),u,d,N_MPC,Q,R,x_ref,ind_sim,Ts,param),[ustar(2:end,:);ustar(end,:)],...
                            [],[],[],[],[],[],@(u)fulltrack_constr_fun(xMPC(:,ind_sim-1),u,d,N,umax,Ts,param),options);
    u                   =   ustar(1,:)';
    [tout,xout]         =   ode45(@(t,x)fulltrackmodel(t,x,u,d,param),[0 Ts],xMPC(:,ind_sim-1));
    xMPC(:,ind_sim)     =   xout(end,:)';
    uMPC(ind_sim-1,:)   =   u;
end

%% Comparison between FHOCP and MPC trajectories
figure(10),subplot(3,1,1),plot(tsim,xMPC(1,:),'b'),grid on, hold on, title('x [km]')
figure(10),subplot(3,1,2),plot(tsim,xMPC(2,:),'b'),grid on, hold on, title('y [km]')
figure(10),subplot(3,1,3),plot(tsim,xMPC(3,:),'b'),grid on, hold on, title('z [km]')

figure(10),subplot(3,1,1),plot(tsim,xsim(1,:),'k'),grid on, hold on, title('x [km]')
figure(10),subplot(3,1,2),plot(tsim,xsim(2,:),'k'),grid on, hold on, title('y [km]')
figure(10),subplot(3,1,3),plot(tsim,xsim(3,:),'k'),grid on, hold on, title('z [km]')

figure(10),subplot(3,1,1),plot(tsim,x_ref(1,1:end-N_MPC),'r--'),grid on, hold on, title('x [km]')
figure(10),subplot(3,1,2),plot(tsim,x_ref(2,1:end-N_MPC),'r--'),grid on, hold on, title('y [km]')
figure(10),subplot(3,1,3),plot(tsim,x_ref(3,1:end-N_MPC),'r--'),grid on, hold on, title('z [km]')

legend('MPC','FHOCP','REF')

figure(11),subplot(3,1,1),plot(tsim,xMPC(7,:),'b'),grid on, hold on, title('alfa-x [rad]')
figure(11),subplot(3,1,2),plot(tsim,xMPC(8,:),'b'),grid on, hold on, title('alfa-y [rad]')
figure(11),subplot(3,1,3),plot(tsim,xMPC(9,:),'b'),grid on, hold on, title('alfa-z [rad]')

figure(11),subplot(3,1,1),plot(tsim,xsim(7,:),'k'),grid on, hold on, title('alfa-x [rad]')
figure(11),subplot(3,1,2),plot(tsim,xsim(8,:),'k'),grid on, hold on, title('alfa-y [rad]')
figure(11),subplot(3,1,3),plot(tsim,xsim(9,:),'k'),grid on, hold on, title('alfa-z [rad]')

figure(11),subplot(3,1,1),plot(tsim,x_ref(7,1:end-N_MPC),'r--'),grid on, hold on, title('alfa-x [rad]')
figure(11),subplot(3,1,2),plot(tsim,x_ref(8,1:end-N_MPC),'r--'),grid on, hold on, title('alfa-y [rad]')
figure(11),subplot(3,1,3),plot(tsim,x_ref(9,1:end-N_MPC),'r--'),grid on, hold on, title('alfa-z [rad]')

legend('MPC','FHOCP','REF')

%% Comparison between FHOCP and MPC control effort
figure(666),subplot(3,1,1),plot(tsim,uMPC(:,1)),grid on, hold on, title('Input acceleration [m/s^2] with MPC vs FHOCP')
figure(666),subplot(3,1,2),plot(tsim,uMPC(:,2)),grid on, hold on
figure(666),subplot(3,1,3),plot(tsim,uMPC(:,3)),grid on, hold on

%% Resulting trajectory with MPC
figure, plot3(xMPC(1,:), xMPC(2,:), xMPC(3,:), 'b', 'LineWidth', 1.5), grid on, hold on
plot3(r_ref(1,1:end-N_MPC),r_ref(2,1:end-N_MPC),r_ref(3,1:end-N_MPC),'r--'),grid on
scatter3(0, 0, 0, 'k', 'filled')
title('MPC fly-around trajectory'), xlabel('x [km]'), ylabel('y [km]'), zlabel('z [km]')
legend('MPC Trajectory', 'Reference Trajectory', 'Target')
axis equal

%% Plot error along the 3-axis
error_x         = xMPC(1,:) - x_ref(1,1:end-N_MPC);
error_y         = xMPC(2,:) - x_ref(2,1:end-N_MPC);
error_z         = xMPC(3,:) - x_ref(3,1:end-N_MPC);
error_anglex    = xMPC(7,:) - x_ref(7,1:end-N_MPC);
error_angley    = xMPC(8,:) - x_ref(8,1:end-N_MPC);
error_anglez    = xMPC(9,:) - x_ref(9,1:end-N_MPC);

e_x         = xsim(1,:) - x_ref(1,1:end-N_MPC);
e_y         = xsim(2,:) - x_ref(2,1:end-N_MPC);
e_z         = xsim(3,:) - x_ref(3,1:end-N_MPC);
e_anglex    = xsim(7,:) - x_ref(7,1:end-N_MPC);
e_angley    = xsim(8,:) - x_ref(8,1:end-N_MPC);
e_anglez    = xsim(9,:) - x_ref(9,1:end-N_MPC);

figure(14),subplot(3,1,1),plot(tsim,error_x,'b','LineWidth', 1.5),grid on,title('Position error (MPC)')
figure(14),subplot(3,1,2),plot(tsim,error_y,'r','LineWidth', 1.5),grid on
figure(14),subplot(3,1,3),plot(tsim,error_z,'g','LineWidth', 1.5),grid on

figure(15),subplot(3,1,1),plot(tsim,error_anglex,'b','LineWidth', 1.5),grid on,title('Angle error (MPC)')
figure(15),subplot(3,1,2),plot(tsim,error_angley,'r','LineWidth', 1.5),grid on
figure(15),subplot(3,1,3),plot(tsim,error_anglez,'g','LineWidth', 1.5),grid on

error_MPC   = [error_x', error_y', error_z'];
error_FHOCP = [e_x', e_y', e_z'];
error_angleMPC   = [error_anglex', error_angley', error_anglez'];
error_angleFHOCP = [e_anglex', e_angley', e_anglez'];
norm_err = zeros(N,4);
for i = 1:N
    norm_err(i,1) = norm(error_MPC(i,:));
    norm_err(i,2) = norm(error_FHOCP(i,:));
    norm_err(i,3) = norm(error_angleMPC(i,:));
    norm_err(i,4) = norm(error_angleFHOCP(i,:));
end

figure(16),subplot(2,1,1),plot(tsim,norm_err(:,1),'r','LineWidth', 1.5),grid on,title('Norm of the position error (MPC)')
figure(16),subplot(2,1,2),plot(tsim,norm_err(:,2),'b','LineWidth', 1.5),grid on,title('Norm of the position error (FHOCP)')
figure(17),subplot(2,1,1),plot(tsim,norm_err(:,3),'r','LineWidth', 1.5),grid on,title('Norm of the angle error (MPC)')
figure(17),subplot(2,1,2),plot(tsim,norm_err(:,4),'b','LineWidth', 1.5),grid on,title('Norm of the angle error (FHOCP)')
