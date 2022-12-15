%% Main file for the MPC control during the tracking of a target position.

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
ic = [-0.01 + 0.02 * randn(1); -0.01 + 0.02 * randn(1); -0.01 + 0.02 * randn(1);
    -0.01*1e-3 + 0.02*1e-3 * randn(1);
    -0.01*1e-3 + 0.02*1e-3 * randn(1); -0.01*1e-3 + 0.02*1e-3 * randn(1)] ;

%% Open loop simulation
Nsim            =   100 ;
x0              =   ic ;
u               =   [1e-3; 0; 0] ;
d               =   [0; 0; 0] ;

xsim            =   zeros(6, Nsim) ;
xsim(:,1)       =   x0 ;
tsim            =   0:Ts:(Nsim-1)*Ts ;

for ind_sim=2:Nsim
    [tout, xout]     = ode45(@(t,x)model(t,x,u,d,param), [0 Ts], xsim(:, ind_sim-1)) ;
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

%% ---------- Reference trajectory (or point) ---------- %
N_MPC       =   40;
N           =   Nsim;
x_ref       =   zeros(6,N+N_MPC);

w_target    =   [0; 0; 0.05];

phi         =   linspace(0,w_target(3)*(N+N_MPC),N+N_MPC);
%r_mag       =   linspace(0.01, 0.003, N);
r_mag0      =   0.003;
r_ref       =   zeros(3,N);
v_ref       =   zeros(3,N);
for i = 1:N+N_MPC
    r_ref(:,i) = [r_mag0 * cos(phi(i)),r_mag0 * sin(phi(i)), 0];
    v_ref(:,i) = cross(w_target,r_ref(:,i));
    
    x_ref(:,i) = [r_ref(:,i); v_ref(:,i)];
end

% x_ref = [0.003; 0; 0; 0; 0; 0];

figure(3),subplot(2,1,1),plot3(x_ref(1,:),x_ref(2,:),x_ref(3,:)),grid on,hold on,title('Reference trajectory')
figure(3),subplot(2,1,1),scatter3(0,0,0,'filled')
figure(3),subplot(2,1,2),plot3(x_ref(4,:),x_ref(5,:),x_ref(6,:)),grid on,title('Reference velocity')

% ----------------------------------------------------- %
close all

%% LQR Control Problem
Q           =   1 ;
R           =   1e-5 ;

K           =   lqr(sys_d,Q,R);

xsim        =   zeros(6,Nsim);
xsim(:,1)   =   x0;
tsim        =   0:Ts:(Nsim-1)*Ts;

ulqr        =   zeros(Nsim,3);
uMax        =   0.02;

for ind_sim=2:Nsim
    u                   = - K * (xsim(:,ind_sim-1)-x_ref(:,ind_sim-1));
    u(u>uMax)           = uMax;
    d                   = [normrnd(0,1e-3); normrnd(0,1e-3); normrnd(0,1e-3)];
    [tout,xout]         = ode45(@(t,x)model(t,x,u,d,param),[0 Ts],xsim(:,ind_sim-1));
    xsim(:,ind_sim)     = xout(end,:)';
    ulqr(ind_sim-1,:)   = u;
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

% close all

figure(6),subplot(2,1,1),plot(tsim,ulqr(:,1)),grid on, hold on, title('Input acceleration [km/s^2] with LQR vs MPC')
figure(6),subplot(2,1,1),plot(tsim,ulqr(:,2)),grid on, hold on
figure(6),subplot(2,1,1),plot(tsim,ulqr(:,3)),grid on, hold on

%% MPC (Close loop simulation w. receding horizon)
Nsim            =   N;
x0              =   ic ;

xMPC            =   zeros(6,Nsim);
uMPC            =   zeros(Nsim,3);
xMPC(:,1)       =   x0;
tsim            =   0:Ts:(Nsim-1)*Ts;
ustar           =   [zeros(N,3);zeros(1,3)];

umax        =   [1e-2 1e-2 1e-2] ;

options     =   optimset('Display','Iter','MaxFunEvals',1e4,'Algorithm','active-set'); 

f = waitbar(0,'1','Name','Running MPC...',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(f,'canceling',0);

for ind_sim=2:Nsim
    % Check for clicked Cancel button
    if getappdata(f,'canceling')
        break
    end
    % Update waitbar and message
    waitbar(ind_sim/Nsim,f,sprintf('%12.9f',ind_sim))

    ustar               =   fmincon(@(u)cw_cost_fun(xMPC(:,ind_sim-1),u,N_MPC,Q,R,Ad,Bd,Cd,Dd,x_ref,ind_sim),[ustar(2:end,:);ustar(end,:)],...
                            [],[],[],[],[],[],@(u)cw_constr_fun(xMPC(:,ind_sim-1),u,N_MPC,Ad,Bd,umax),options);
    u                   =   ustar(1,:)';
    d                   =   [normrnd(0,1e-3); normrnd(0,1e-3); normrnd(0,1e-3)];
    [tout,xout]         =   ode45(@(t,x)model(t,x,u,d,param),[0 Ts],xMPC(:,ind_sim-1));
    xMPC(:,ind_sim)     =   xout(end,:)';
    uMPC(ind_sim-1,:)   =   u;
end

delete(f)

figure(1),subplot(4,1,1),plot(tsim,xMPC(1,:)),grid on, hold on, title('x [km]')
figure(1),subplot(4,1,2),plot(tsim,xMPC(2,:)),grid on, hold on, title('y [km]')
figure(1),subplot(4,1,3),plot(tsim,xMPC(3,:)),grid on, hold on, title('z [km]')

%% Comparison between LQR and MPC trajectories
figure(7),subplot(3,1,1),plot(tsim,xMPC(1,:),'b'),grid on, hold on, title('x [km]')
figure(7),subplot(3,1,2),plot(tsim,xMPC(2,:),'b'),grid on, hold on, title('y [km]')
figure(7),subplot(3,1,3),plot(tsim,xMPC(3,:),'b'),grid on, hold on, title('z [km]')

figure(7),subplot(3,1,1),plot(tsim,xsim(1,:),'k'),grid on, hold on, title('x [km]')
figure(7),subplot(3,1,2),plot(tsim,xsim(2,:),'k'),grid on, hold on, title('y [km]')
figure(7),subplot(3,1,3),plot(tsim,xsim(3,:),'k'),grid on, hold on, title('z [km]')

figure(7),subplot(3,1,1),plot(x_ref(1,1:end-N_MPC),'r--'),grid on, hold on, title('x [km]')
figure(7),subplot(3,1,2),plot(x_ref(2,1:end-N_MPC),'r--'),grid on, hold on, title('y [km]')
figure(7),subplot(3,1,3),plot(x_ref(3,1:end-N_MPC),'r--'),grid on, hold on, title('z [km]')

legend('MPC', 'LQR', 'REF')

%% Comparison between FHOCP and MPC control effort
figure(6),subplot(2,1,2),plot(tsim,uMPC(:,1)),grid on, hold on, title('Input acceleration [km/s^2] with LQR vs MPC')
figure(6),subplot(2,1,2),plot(tsim,uMPC(:,2)),grid on, hold on
figure(6),subplot(2,1,2),plot(tsim,uMPC(:,3)),grid on, hold on

legend('LQR', 'MPC')

%% Resulting trajectory with MPC
figure, plot3(xMPC(1,:), xMPC(2,:), xMPC(3,:), 'b', 'LineWidth', 1.5), grid on, hold on
plot3(x_ref(1,:),x_ref(2,:),x_ref(3,:),'r--'),grid on
scatter3(0, 0, 0, 'k', 'filled')
title('MPC fly-around trajectory'), xlabel('x [km]'), ylabel('y [km]'), zlabel('z [km]')
legend('MPC Trajectory', 'Reference Trajectory', 'Target')
axis equal

%% Resulting trajectory with LQR
figure, plot3(xsim(1,:), xsim(2,:), xsim(3,:), 'k', 'LineWidth', 1.5), grid on, hold on
plot3(x_ref(1,1:end-N_MPC),x_ref(2,1:end-N_MPC),x_ref(3,1:end-N_MPC),'r--'),grid on
scatter3(0, 0, 0, 'k', 'filled')
title('LQR fly-around trajectory'), xlabel('x [km]'), ylabel('y [km]'), zlabel('z [km]')
legend('LQR Trajectory', 'Reference Trajectory', 'Target')
% axis equal

%% Plot error along the 3-axis
error_x = xMPC(1,:) - x_ref(1,1:end-N_MPC);
error_y = xMPC(2,:) - x_ref(2,1:end-N_MPC);
error_z = xMPC(3,:) - x_ref(3,1:end-N_MPC);

figure,subplot(3,1,1),plot(tsim,error_x,'b','LineWidth', 1.5),grid on
subplot(3,1,2),plot(tsim,error_y,'r','LineWidth', 1.5),grid on
subplot(3,1,3),plot(tsim,error_z,'g','LineWidth', 1.5),grid on

error = [error_x', error_y', error_z'];
norm_err = zeros(N,1);
for i = 1:N
    norm_err(i) = norm(error(i,:));
end

figure,plot(tsim,norm_err,'r','LineWidth', 1.5),grid on,title('Norm of the position error')
