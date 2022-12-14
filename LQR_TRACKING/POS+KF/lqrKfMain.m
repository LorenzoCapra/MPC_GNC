%% Main file for the LQR control during the tracking of a target position using Kalman Filter to estimate the state.

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
Ts = .1 ; % [s]

% Initial conditions:
ic = [-0.1 + 0.2 * randn(1); -0.1 + 0.2 * randn(1); -0.1 + 0.2 * randn(1); ...
    -0.01*1e-3 + 0.02*1e-3 * randn(1); -0.01*1e-3 + 0.02*1e-3 * randn(1); -0.01*1e-3 + 0.02*1e-3 * randn(1)] ;

%% State space model + Kalman Filter
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

C = eye(6);

D = zeros(size(C,1),size(Bu,2)) ;

system          =   ss(A,Bu,C,D);

sys_d           =   c2d(system,Ts,'zoh');
[Ad,Bd,Cd,Dd]   =   ssdata(sys_d);

%% Open loop simulation for the nominal system
Nsim            =   1000 ;
t               =   0:Ts:Nsim*Ts;
x0              =   ic ;
u0              =   [1e-3*sin(t); -1e-3*sin(t); 5e-4*sin(t)] ;

% Simulate with control input
[y,ti]           =   lsim(system,u0,t,x0);
figure,plot(ti,y(:,1:3),'Linewidth',1.5),grid on,hold on

[yd,ti]           =   lsim(sys_d,u0,t,x0);
plot(ti,yd(:,1:3),'--'),grid on,hold on

%% Open loop simulation for the nominal + noise system
Vd = .1*eye(6);   % disturbance covariance
Vn = .2*eye(size(C,1));   % noise covariance
param.Vd = Vd;
param.Vn = Vn;

Bf = [Bu, Vd, 0*Vn];

sysNoise = ss(A,Bf,C,[zeros(6,9),Vn]);

uDist = randn(6,size(t,2));
uNoise = randn(6,size(t,2));
u0 = [1e-3* sin(t); -1e-3* sin(t); 1e-3* sin(t)] ;

uAug = [u0; Vd*Vd*uDist; Vn*Vn*uNoise];

[y_noise,ti,xn] = lsim(sysNoise, uAug, t,x0);
figure,plot(ti,y_noise(:,1:3)),grid on, hold on
% plot(ti,xn(:,1:3)),grid on,hold on

%% Open loop simulation for the nominal system + disturbance
sysDist = ss(A,Bf,eye(6),0);

[ydist,ti] = lsim(sysDist,uAug,t,x0);
plot(ti,ydist(:,1:3),'k--'),grid on, hold on

%% Build Kalman Filter
[L,Pe,E] = lqe(A,Vd,C,Vd,Vn); % design kalman filter

sysKf = ss(A-L*C, [Bd,L], C, 0*[Bu,L]); % Kalman filter estimator

[yKF,ti] = lsim(sysKf,[u0;y_noise'],t,x0);
plot(ti,yKF(:,1:3),'g','LineWidth',2)

xlabel('Time [s]'), ylabel('[km]'), title('Open loop simulation with KF')
legend('Noisy Measurement', '','','Model Measurment', '','','Filtered Measurement')

%% Apply Kalman Filter during the simulation
xsim            =   zeros(size(ic,1),Nsim);
xsim(:,1)       =   x0;
xdsim           =   zeros(size(ic,1),Nsim);
xdsim(:,1)      =   x0;
ysim            =   zeros(size(ic,1),Nsim);
ysim(:,1)       =   C*x0;
xHatsim         =   zeros(size(ic,1),Nsim);
xHatsim(:,1)    =   x0;
P               =   eye(size(ic,1));
tsim            =   0:Ts:(Nsim-1)*Ts;

for ind_sim = 2:Nsim
    u               = randn(3,1)*1e-4;
    d               = [0;0;0];
    % Integrate nominal dynamics
    [tout1, xout]   = ode45(@(t,x)model(t,x,u,d,param), [0 Ts], xsim(:, ind_sim-1)) ;
    xsim(:,ind_sim) = xout(end,:)' ;
    ysim(:,ind_sim) = system.C * xsim(:,ind_sim);

    % Integrate disturbed dynamic measurements
    uAug = [u; Vd*Vd*randn(6,1); Vn*Vn*randn(6,1)];

    [tout2, xdout]   = ode45(@(t,x)modelDist(t,x,uA,sysNoise), [0 Ts], xdsim(:, ind_sim-1)) ;
    xdsim(:,ind_sim) = sysNoise.C*xdout(end,:)';

    % Estimate the state
    [tout3, xkout]   = ode45(@(t,x)modelKalman(t,x,u,ysim(:,ind_sim),sysKf), [0 Ts], xHatsim(:, ind_sim-1)) ;
    xHatsim(:,ind_sim)  = xkout(end,:)';
end

%% Plot the state evolution and reconstruction
figure, plot(tsim, xsim(1:3,:), 'k--'),grid on,hold on
% plot(tsim, xsim(1:3,:)+.01*randn(3,Nsim)),grid on, hold on
plot(tsim, xdsim(1:3,:)),grid on, hold on
plot(tsim, xHatsim(1:3,:), 'g', 'LineWidth',1.5),grid on, hold on
xlabel('Time [s]'), ylabel('[km]'), title('Kalman Filter application during simulation')

legend('Model Measurment', '','','Noisy Measurement', '','','Filtered Measurement')

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

xsim            =   zeros(6,Nsim);
xsim(:,1)       =   x0;
xHatsim         =   zeros(6,Nsim);
xHatsim(:,1)    =   x0;
tsim            =   0:Ts:(Nsim-1)*Ts;

ulqr        =   zeros(Nsim,3);
uMax        =   0.02;

for ind_sim=2:Nsim
    u                   = - K * (xHatsim(:,ind_sim-1)-x_ref(:,ind_sim-1));
    u(u>uMax)           = uMax;
%     d                   = [normrnd(0,1e-3); normrnd(0,1e-3); normrnd(0,1e-3)];
%     [tout,xout]         = ode45(@(t,x)modelDist(t,x,u,param,sysFullOutput),[0 Ts],xsim(:,ind_sim-1));
    [xout,ti] = lsim(sysFullOutput,[u,u; Vd*Vd*randn(6,1),Vd*Vd*randn(6,1); Vn*Vn*randn(6,1),Vn*Vn*randn(6,1)],[0,Ts]);
    xsim(:,ind_sim)     = xout(end,:)';
    ulqr(ind_sim-1,:)   = u;

    xHatsim(:,ind_sim)  = (sysKf.A * xHatsim(:,ind_sim-1) + sysKf.B * [u; xsim(:,ind_sim)]) * Ts ;
end

% Position vs reference
figure(4),subplot(3,1,1),plot(tsim,xsim(1,:)),grid on, hold on, title('x [km]')
figure(4),subplot(3,1,2),plot(tsim,xsim(2,:)),grid on, hold on, title('y [km]')
figure(4),subplot(3,1,3),plot(tsim,xsim(3,:)),grid on, hold on, title('z [km]')

figure(4),subplot(3,1,1),plot(tsim,xHatsim(1,:)),grid on, hold on, title('x [km]')
figure(4),subplot(3,1,2),plot(tsim,xHatsim(2,:)),grid on, hold on, title('y [km]')
figure(4),subplot(3,1,3),plot(tsim,xHatsim(3,:)),grid on, hold on, title('z [km]')

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

%% %% MPC (Close loop simulation w. receding horizon)
% Nsim            =   N;
% x0              =   ic ;
% 
% xMPC            =   zeros(6,Nsim);
% uMPC            =   zeros(Nsim,3);
% xMPC(:,1)       =   x0;
% tsim            =   0:Ts:(Nsim-1)*Ts;
% ustar           =   [zeros(N,3);zeros(1,3)];
% 
% umax        =   [1e-2 1e-2 1e-2] ;
% 
% options     =   optimset('Display','Iter','MaxFunEvals',1e4,'Algorithm','active-set'); 
% 
% f = waitbar(0,'1','Name','Running MPC...',...
%     'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
% setappdata(f,'canceling',0);
% 
% for ind_sim=2:Nsim
%     % Check for clicked Cancel button
%     if getappdata(f,'canceling')
%         break
%     end
%     % Update waitbar and message
%     waitbar(ind_sim/Nsim,f,sprintf('%12.9f',ind_sim))
% 
%     ustar               =   fmincon(@(u)cw_cost_fun(xMPC(:,ind_sim-1),u,N_MPC,Q,R,Ad,Bd,Cd,Dd,x_ref,ind_sim),[ustar(2:end,:);ustar(end,:)],...
%                             [],[],[],[],[],[],@(u)cw_constr_fun(xMPC(:,ind_sim-1),u,N_MPC,Ad,Bd,umax),options);
%     u                   =   ustar(1,:)';
%     d                   =   [normrnd(0,1e-3); normrnd(0,1e-3); normrnd(0,1e-3)];
%     [tout,xout]         =   ode45(@(t,x)model(t,x,u,d,param),[0 Ts],xMPC(:,ind_sim-1));
%     xMPC(:,ind_sim)     =   xout(end,:)';
%     uMPC(ind_sim-1,:)   =   u;
% end
% 
% delete(f)
% 
% figure(1),subplot(4,1,1),plot(tsim,xMPC(1,:)),grid on, hold on, title('x [km]')
% figure(1),subplot(4,1,2),plot(tsim,xMPC(2,:)),grid on, hold on, title('y [km]')
% figure(1),subplot(4,1,3),plot(tsim,xMPC(3,:)),grid on, hold on, title('z [km]')
% 
% %% Comparison between LQR and MPC trajectories
% figure(7),subplot(3,1,1),plot(tsim,xMPC(1,:),'b'),grid on, hold on, title('x [km]')
% figure(7),subplot(3,1,2),plot(tsim,xMPC(2,:),'b'),grid on, hold on, title('y [km]')
% figure(7),subplot(3,1,3),plot(tsim,xMPC(3,:),'b'),grid on, hold on, title('z [km]')
% 
% figure(7),subplot(3,1,1),plot(tsim,xsim(1,:),'k'),grid on, hold on, title('x [km]')
% figure(7),subplot(3,1,2),plot(tsim,xsim(2,:),'k'),grid on, hold on, title('y [km]')
% figure(7),subplot(3,1,3),plot(tsim,xsim(3,:),'k'),grid on, hold on, title('z [km]')
% 
% figure(7),subplot(3,1,1),plot(x_ref(1,1:end-N_MPC),'r--'),grid on, hold on, title('x [km]')
% figure(7),subplot(3,1,2),plot(x_ref(2,1:end-N_MPC),'r--'),grid on, hold on, title('y [km]')
% figure(7),subplot(3,1,3),plot(x_ref(3,1:end-N_MPC),'r--'),grid on, hold on, title('z [km]')
% 
% legend('MPC', 'LQR', 'REF')
% 
% %% Comparison between FHOCP and MPC control effort
% figure(6),subplot(2,1,2),plot(tsim,uMPC(:,1)),grid on, hold on, title('Input acceleration [km/s^2] with LQR vs MPC')
% figure(6),subplot(2,1,2),plot(tsim,uMPC(:,2)),grid on, hold on
% figure(6),subplot(2,1,2),plot(tsim,uMPC(:,3)),grid on, hold on
% 
% legend('LQR', 'MPC')
% 
% %% Resulting trajectory with MPC
% figure, plot3(xMPC(1,:), xMPC(2,:), xMPC(3,:), 'b', 'LineWidth', 1.5), grid on, hold on
% plot3(x_ref(1,:),x_ref(2,:),x_ref(3,:),'r--'),grid on
% scatter3(0, 0, 0, 'k', 'filled')
% title('MPC fly-around trajectory'), xlabel('x [km]'), ylabel('y [km]'), zlabel('z [km]')
% legend('MPC Trajectory', 'Reference Trajectory', 'Target')
% axis equal
% 
% %% Resulting trajectory with LQR
% figure, plot3(xsim(1,:), xsim(2,:), xsim(3,:), 'k', 'LineWidth', 1.5), grid on, hold on
% plot3(x_ref(1,1:end-N_MPC),x_ref(2,1:end-N_MPC),x_ref(3,1:end-N_MPC),'r--'),grid on
% scatter3(0, 0, 0, 'k', 'filled')
% title('LQR fly-around trajectory'), xlabel('x [km]'), ylabel('y [km]'), zlabel('z [km]')
% legend('LQR Trajectory', 'Reference Trajectory', 'Target')
% % axis equal
% 
% %% Plot error along the 3-axis
% error_x = xMPC(1,:) - x_ref(1,1:end-N_MPC);
% error_y = xMPC(2,:) - x_ref(2,1:end-N_MPC);
% error_z = xMPC(3,:) - x_ref(3,1:end-N_MPC);
% 
% figure,subplot(3,1,1),plot(tsim,error_x,'b','LineWidth', 1.5),grid on
% subplot(3,1,2),plot(tsim,error_y,'r','LineWidth', 1.5),grid on
% subplot(3,1,3),plot(tsim,error_z,'g','LineWidth', 1.5),grid on
% 
% error = [error_x', error_y', error_z'];
% norm_err = zeros(N,1);
% for i = 1:N
%     norm_err(i) = norm(error(i,:));
% end
% 
% figure,plot(tsim,norm_err,'r','LineWidth', 1.5),grid on,title('Norm of the position error')
