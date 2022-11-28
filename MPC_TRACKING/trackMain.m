%% Main file for the MPC control during the tracking of a target pose.

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
Ts = 0.5 ; % [s]

% Initial conditions:
w_target    =   [0.0; 0.0; 0.05];
r_mag0      =   0.003;
r0          =   [r_mag0; 0; 0];
v0          =   cross(w_target, r0);
ic          =   [r0; v0] ;
%ic = [0.01; 0; 0; v0];

%% Open loop simulation
Nsim            =   100 ;
x0              =   ic ;
u               =   [1e-3; 0; 0] ;
d               =   [normrnd(0, 1e-4); normrnd(0, 1e-4); normrnd(0, 1e-4)];

xsim            =   zeros(6, Nsim) ;
xsim(:,1)       =   x0 ;
tsim            =   0:Ts:(Nsim-1)*Ts ;

for ind_sim=2:Nsim
    [tout, xout]     = ode45(@(t,x)trackmodel(t,x,u,d,param), [0 Ts], xsim(:, ind_sim-1)) ;
    xsim(:,ind_sim)  = xout(end,:)' ;
end

figure(1),subplot(4,1,1),plot(tsim,xsim(1,:)),grid on, hold on, title('x [km]')
figure(1),subplot(4,1,2),plot(tsim,xsim(2,:)),grid on, hold on, title('y [km]')
figure(1),subplot(4,1,3),plot(tsim,xsim(3,:)),grid on, hold on, title('z [km]')
figure(1),subplot(4,1,4),plot(tsim,u*ones(1,Nsim)),grid on, hold on, title('Input acceleration [km/s^2]')

figure(2),plot3(xsim(1,:),xsim(2,:),xsim(3,:),'LineWidth',1.5),grid on,hold on,
title('Relative Motion Chaser-Target')
xlabel('x [km]'), ylabel('y [km]'), zlabel('z [km]')
figure(2),scatter3(0,0,0, 'filled'),grid on

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
x0          =   ic ;
N           =   120 ;
Q           =   1 ;
R           =   1e-5 ;
umax        =   [1e-2 1e-2 1e-2] ;

% ---------- Reference trajectory (or point) ---------- %
N_MPC       =   10;
x_ref       =   zeros(6,N+N_MPC);

phi         =   linspace(0,w_target(3)*N,N+N_MPC);
%r_mag       =   linspace(0.01, 0.003, N);
r_ref       =   zeros(3,N+N_MPC);
v_ref       =   zeros(3,N+N_MPC);
for i = 1:N+N_MPC
    r_ref(:,i) = [r_mag0 * cos(phi(i)),r_mag0 * sin(phi(i)), 0];
    v_ref(:,i) = cross(w_target,r_ref(:,i));
    
    x_ref(:,i) = [r_ref(:,i); v_ref(:,i)];
end

figure(3),subplot(2,1,1),plot3(r_ref(1,:),r_ref(2,:),r_ref(3,:)),grid on,hold on,title('Reference trajectory')
figure(3),subplot(2,1,1),scatter3(0,0,0,'filled')
figure(3),subplot(2,1,2),plot3(v_ref(1,:),v_ref(2,:),v_ref(3,:)),grid on,title('Reference velocity')

% ----------------------------------------------------- %

%% Set the options and solve the optimization problem in open loop
options     =   optimset('Display','Iter','MaxFunEvals',1e4,'Algorithm','active-set');  

ustar   = fmincon(@(u)track_cost_fun(x0,u,N,Q,R,Ad,Bd,Cd,Dd,x_ref,1),zeros(N,3),...
    [],[],[],[],[],[],@(u)track_constr_fun(x0,u,N,Ad,Bd,umax),options);

[Constr,Constr_eq,xpred]=track_constr_fun(x0,ustar,N,Ad,Bd,umax);

%% Open loop simulation - FHOCP solution
Nsim            =   N;
d               =   [0;0;0];

xsim            =   zeros(6,Nsim);
xsim(:,1)       =   x0;
tsim            =   0:Ts:(Nsim-1)*Ts;

for ind_sim=2:Nsim
    u               = ustar(ind_sim-1,:)';
    [tout,xout]     = ode45(@(t,x)trackmodel(t,x,u,d,param),[0 Ts],xsim(:,ind_sim-1));
    xsim(:,ind_sim) = xout(end,:)';
end

% Position vs reference
figure(4),subplot(3,1,1),plot(tsim,xsim(1,:)),grid on, hold on, title('x [km]')
figure(4),subplot(3,1,2),plot(tsim,xsim(2,:)),grid on, hold on, title('y [km]')
figure(4),subplot(3,1,3),plot(tsim,xsim(3,:)),grid on, hold on, title('z [km]')

figure(4),subplot(3,1,1),plot(tsim,x_ref(1,1:end-N_MPC),'r--'),grid on, hold on
figure(4),subplot(3,1,2),plot(tsim,x_ref(2,1:end-N_MPC),'r--'),grid on, hold on
figure(4),subplot(3,1,3),plot(tsim,x_ref(3,1:end-N_MPC),'r--'),grid on, hold on

% Velocity vs reference
figure(5),subplot(3,1,1),plot(tsim,xsim(4,:)),grid on, hold on, title('x [km/s]')
figure(5),subplot(3,1,2),plot(tsim,xsim(5,:)),grid on, hold on, title('y [km/s]')
figure(5),subplot(3,1,3),plot(tsim,xsim(6,:)),grid on, hold on, title('z [km/s]')

figure(5),subplot(3,1,1),plot(tsim,x_ref(4,1:end-N_MPC),'r--'),grid on, hold on
figure(5),subplot(3,1,2),plot(tsim,x_ref(5,1:end-N_MPC),'r--'),grid on, hold on
figure(5),subplot(3,1,3),plot(tsim,x_ref(6,1:end-N_MPC),'r--'),grid on, hold on

% Plot control effort
figure(6),subplot(3,1,1),plot(tsim,ustar(:,1)),grid on, hold on,
figure(6),subplot(3,1,2),plot(tsim,ustar(:,2)),grid on, hold on,
figure(6),subplot(3,1,3),plot(tsim,ustar(:,3)),grid on, hold on

%% MPC (Close loop simulation w. receding horizon)
Nsim            =   N;
% N_MPC           =   10;
x0              =   ic ;
d               =   [0;0;0];

xMPC            =   zeros(6,Nsim);
uMPC            =   zeros(Nsim,3);
xMPC(:,1)       =   x0;
tsim            =   0:Ts:(Nsim-1)*Ts;
ustar           =   [zeros(N,3);zeros(1,3)];

for ind_sim=2:Nsim
    ustar               =   fmincon(@(u)track_cost_fun(xMPC(:,ind_sim-1),u,N_MPC,Q,R,Ad,Bd,Cd,Dd,x_ref,ind_sim),[ustar(2:end,:);ustar(end,:)],...
                            [],[],[],[],[],[],@(u)track_constr_fun(xMPC(:,ind_sim-1),u,N_MPC,Ad,Bd,umax),options);
    u                   =   ustar(1,:)';
    [tout,xout]         =   ode45(@(t,x)trackmodel(t,x,u,d,param),[0 Ts],xMPC(:,ind_sim-1));
    xMPC(:,ind_sim)     =   xout(end,:)';
    uMPC(ind_sim-1,:)   =   u;
end

%% Comparison between FHOCP and MPC trajectories
figure(7),subplot(3,1,1),plot(tsim,xMPC(1,:),'b'),grid on, hold on, title('x [km]')
figure(7),subplot(3,1,2),plot(tsim,xMPC(2,:),'b'),grid on, hold on, title('y [km]')
figure(7),subplot(3,1,3),plot(tsim,xMPC(3,:),'b'),grid on, hold on, title('z [km]')

figure(7),subplot(3,1,1),plot(tsim,xsim(1,:),'k'),grid on, hold on, title('x [km]')
figure(7),subplot(3,1,2),plot(tsim,xsim(2,:),'k'),grid on, hold on, title('y [km]')
figure(7),subplot(3,1,3),plot(tsim,xsim(3,:),'k'),grid on, hold on, title('z [km]')

figure(7),subplot(3,1,1),plot(tsim,x_ref(1,1:end-N_MPC),'r--'),grid on, hold on, title('x [km]')
figure(7),subplot(3,1,2),plot(tsim,x_ref(2,1:end-N_MPC),'r--'),grid on, hold on, title('y [km]')
figure(7),subplot(3,1,3),plot(tsim,x_ref(3,1:end-N_MPC),'r--'),grid on, hold on, title('z [km]')

legend('MPC','FHOCP','REF')

%% Comparison between FHOCP and MPC control effort
figure(6),subplot(3,1,1),plot(tsim,uMPC(:,1)),grid on, hold on, title('Input acceleration [m/s^2] with MPC vs FHOCP')
figure(6),subplot(3,1,2),plot(tsim,uMPC(:,2)),grid on, hold on
figure(6),subplot(3,1,3),plot(tsim,uMPC(:,3)),grid on, hold on

%% Resulting trajectory with MPC
figure, plot3(xMPC(1,:), xMPC(2,:), xMPC(3,:), 'b', 'LineWidth', 1.5), grid on, hold on
plot3(r_ref(1,1:end-N_MPC),r_ref(2,1:end-N_MPC),r_ref(3,1:end-N_MPC),'r--'),grid on
scatter3(0, 0, 0, 'k', 'filled')
title('MPC fly-around trajectory'), xlabel('x [km]'), ylabel('y [km]'), zlabel('z [km]')
legend('MPC Trajectory', 'Reference Trajectory', 'Target')
axis equal

%% Plot error along the 3-axis
error_x = xMPC(1,:) - x_ref(1,1:end-N_MPC);
error_y = xMPC(2,:) - x_ref(2,1:end-N_MPC);
error_z = xMPC(3,:) - x_ref(3,1:end-N_MPC);

e_x = xsim(1,:) - x_ref(1,1:end-N_MPC);
e_y = xsim(2,:) - x_ref(2,1:end-N_MPC);
e_z = xsim(3,:) - x_ref(3,1:end-N_MPC);

figure,subplot(3,1,1),plot(tsim,error_x,'b','LineWidth', 1.5),grid on
subplot(3,1,2),plot(tsim,error_y,'r','LineWidth', 1.5),grid on
subplot(3,1,3),plot(tsim,error_z,'g','LineWidth', 1.5),grid on

error_MPC   = [error_x', error_y', error_z'];
error_FHOCP = [e_x', e_y', e_z'];
norm_err = zeros(N,2);
for i = 1:N
    norm_err(i,1) = norm(error_MPC(i,:));
    norm_err(i,2) = norm(error_FHOCP(i,:));
end

figure(9),subplot(2,1,1),plot(tsim,norm_err(:,1),'r','LineWidth', 1.5),grid on,title('Norm of the position error (MPC)')
figure(9),subplot(2,1,2),plot(tsim,norm_err(:,2),'b','LineWidth', 1.5),grid on,title('Norm of the position error (FHOCP)')