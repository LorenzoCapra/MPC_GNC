%% Main file for the MPC control of the chaser trajectory around a target.

clc
clear
close all

%% System parameters
param.mu = 398600.0 ;                   % Earth gravitational parameter
param.Re = 6378.0 ;                     % Earth radius
param.h = 500.0 ;                       % Orbit altitude
param.R = param.Re + param.h ;          % Orbit radius
param.n = sqrt(param.mu/(param.R^3)) ;  % Mean orbital motion

% Sampling time
Ts = 1 ; % [s]

% Initial conditions:
ic = [-0.1 + 0.2 * randn(1); -0.1 + 0.2 * randn(1); -0.1 + 0.2 * randn(1); ...
    -0.01*1e-3 + 0.02*1e-3 * randn(1); -0.01*1e-3 + 0.02*1e-3 * randn(1); -0.01*1e-3 + 0.02*1e-3 * randn(1)] ;

%% Open loop simulation
Nsim            =   50 ;
x0              =   zeros(6, 1) + ic ; %[0.01; -0.1; 0.2; 0.1*1e-3; -0.02*1e-3; 0.001*1e-3] ;
u               =   [1e-3; 0; 0] ;
d               =   [0; 0; 0] ;

xsim            =   zeros(6, Nsim) ;
xsim(:,1)       =   x0 ;
tsim            =   0:Ts:(Nsim-1)*Ts ;

for ind_sim=2:Nsim
    [tout, xout]     = ode45(@(t,x)cwmodel(t,x,u,d,param), [0 Ts], xsim(:, ind_sim-1)) ;
    xsim(:,ind_sim) = xout(end,:)' ;
end

figure(1),subplot(4,1,1),plot(tsim,xsim(1,:)),grid on, hold on, title('X [km]')
figure(1),subplot(4,1,2),plot(tsim,xsim(2,:)),grid on, hold on, title('Y [km]')
figure(1),subplot(4,1,3),plot(tsim,xsim(3,:)),grid on, hold on, title('Z [km]')
figure(1),subplot(4,1,4),plot(tsim,u*ones(1,Nsim)),grid on, hold on, title('Input acceleration [m/s^2]')

%% Discrete time model
n = param.n ;

A = [0, 0, 0, 1, 0, 0;
     0, 0, 0, 0, 1, 0;
     0, 0, 0, 0, 0, 1;
     3*n^2, 0, 0, 0, 2*n, 0;
     0, 0, 0, -2*n, 0, 0;
     0, 0, -n^2, 0, 0, 0] ;
 
Bu = [0, 0, 0; 
      0, 0, 0;
      0, 0, 0;
      1, 0, 0;
      0, 1, 0;
      0, 0, 1] ;

C = [1, 0, 0, 0, 0, 0;
     0, 1, 0, 0, 0, 0;
     0, 0, 1, 0, 0, 0;
     0, 0, 0, 1, 0, 0;
     0, 0, 0, 0, 1, 0;
     0, 0, 0, 0, 0, 1];

D = 0 ;

system          =   ss(A,Bu,C,D);
sys_d           =   c2d(system,Ts,'tustin');
[Ad,Bd,Cd,Dd]   =   ssdata(sys_d);

%% Finite Horizon Optimal Control Problem
x0          =   zeros(6, 1) + ic ; %[0.01; -0.1; 0.2; 0.1*1e-3; -0.02*1e-3; 0.001*1e-3] ;
N           =   60 ;
Q           =   1 ;
R           =   1e-5 ;
umax        =   [1e-2 1e-2 1e-2] ;
x_ref       =   [0.0; 0.0; -0.003; 0.0; 0.0; 0.0] ;

options     =   optimset('Display','Iter','MaxFunEvals',1e4,'Algorithm','active-set');  

ustar   = fmincon(@(u)cw_cost_fun(x0,u,N,Q,R,Ad,Bd,Cd,Dd,x_ref),zeros(N,3),...
    [],[],[],[],[],[],@(u)cw_constr_fun(x0,u,N,Ad,Bd,umax),options);

figure,plot(ustar)

[Constr,Constr_eq,xpred]=cw_constr_fun(x0,ustar,N,Ad,Bd,umax);
% figure,plot(0:Ts:N*Ts,xpred(2,:))
% figure,plot(0:Ts:(N-1)*Ts,Torsion)

%% Open loop simulation - FHOCP solution
Nsim            =   N;
d               =   [normrnd(0, 1e-4); normrnd(0, 1e-4); normrnd(0, 1e-4)];

xsim            =   zeros(6,Nsim);
xsim(:,1)       =   x0;
tsim            =   0:Ts:(Nsim-1)*Ts;

for ind_sim=2:Nsim
    u               = ustar(ind_sim-1,:)';
    [tout,xout]     = ode45(@(t,x)cwmodel(t,x,u,d,param),[0 Ts],xsim(:,ind_sim-1));
    xsim(:,ind_sim) = xout(end,:)';
end

figure(1),subplot(4,1,1),plot(tsim,xsim(1,:)),grid on, hold on, title('X [km]')
figure(1),subplot(4,1,2),plot(tsim,xsim(2,:)),grid on, hold on, title('Y [km]')
figure(1),subplot(4,1,3),plot(tsim,xsim(3,:)),grid on, hold on, title('Z [km]')
figure(1),subplot(4,1,4),plot(tsim,ustar*ones(3,Nsim)),grid on, hold on, title('Input acceleration [m/s^2]')

figure(3),subplot(2,1,1),plot(tsim,ustar*ones(3,Nsim)),grid on, hold on, title('Input acceleration [m/s^2] with FHOCP')

%% MPC (Close loop simulation w. receding horizon)
Nsim            =   N;
N_MPC           =   50;
x0              =   zeros(6,1) + ic ; %[0.01; -0.1; 0.2; 0.1*1e-3; -0.02*1e-3; 0.001*1e-3];
d               =   [1e-4; -2.3e-5; 7e-6];

xMPC            =   zeros(6,Nsim);
uMPC            =   zeros(Nsim,3);
xMPC(:,1)       =   x0;
tsim            =   0:Ts:(Nsim-1)*Ts;
ustar           =   [zeros(N,3);zeros(1,3)];

for ind_sim=2:Nsim
    ustar               =   fmincon(@(u)cw_cost_fun(xMPC(:,ind_sim-1),u,N_MPC,Q,R,Ad,Bd,Cd,Dd,x_ref),[ustar(2:end,:);ustar(end,:)],...
                        [],[],[],[],[],[],@(u)cw_constr_fun(xMPC(:,ind_sim-1),u,N_MPC,Ad,Bd,umax),options);
    u                   =   ustar(1,:)';
    [tout,xout]         =   ode45(@(t,x)cwmodel(t,x,u,d,param),[0 Ts],xMPC(:,ind_sim-1));
    xMPC(:,ind_sim)     =   xout(end,:)';
    uMPC(ind_sim-1,:)   =   u;
end

figure(1),subplot(4,1,1),plot(tsim,xMPC(1,:)),grid on, hold on, title('X [km]')
figure(1),subplot(4,1,2),plot(tsim,xMPC(2,:)),grid on, hold on, title('Y [km]')
figure(1),subplot(4,1,3),plot(tsim,xMPC(3,:)),grid on, hold on, title('Z [km]')
figure(1),subplot(4,1,4),plot(tsim,uMPC*ones(3,Nsim)),grid on, hold on, title('Input acceleration [m/s^2]')

%% Comparison between FHOCP and MPC trajectories
figure(2),subplot(3,1,1),plot(tsim,xMPC(1,:)),grid on, hold on, title('X [km]')
figure(2),subplot(3,1,2),plot(tsim,xMPC(2,:)),grid on, hold on, title('Y [km]')
figure(2),subplot(3,1,3),plot(tsim,xMPC(3,:)),grid on, hold on, title('Z [km]')

figure(2),subplot(3,1,1),plot(tsim,xsim(1,:)),grid on, hold on, title('X [km]')
figure(2),subplot(3,1,2),plot(tsim,xsim(2,:)),grid on, hold on, title('Y [km]')
figure(2),subplot(3,1,3),plot(tsim,xsim(3,:)),grid on, hold on, title('Z [km]')

%% Comparison between FHOCP and MPC control effort
figure(3),subplot(2,1,2),plot(tsim,uMPC*ones(3,Nsim)),grid on, hold on, title('Input acceleration [m/s^2] with MPC')

%% Resulting trajectory with MPC
figure, plot3(xMPC(1,:), xMPC(2,:), xMPC(3,:), 'r', 'LineWidth', 1.5), grid on, hold on
scatter3(0, 0, 0, 'k', 'filled')
title('MPC rendezvous trajectory'), xlabel('x [km]'), ylabel('y [km]'), zlabel('z [km]')
legend('Trajectory', 'Target')
