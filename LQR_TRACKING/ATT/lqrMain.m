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
w_chaser    =   [0; 0; 0];

ic          =   [w_chaser;o0] ;

%% Open loop simulation
Nsim        =   100 ;

x0          =   ic ;

u           =   [normrnd(0, 1e-3); normrnd(0, 1e-3); normrnd(0, 1e-3)] ;
d           =   [normrnd(0, 1e-4); normrnd(0, 1e-4); normrnd(0, 1e-4)] ;

xsim        =   zeros(size(ic, 1), Nsim) ;
xsim(:,1)   =   x0 ;

tsim        =   0:Ts:(Nsim-1)*Ts ;

for ind_sim=2:Nsim
    xn               = linearAdcsModel(xsim(:,ind_sim-1),u,d,param);
    xsim(:,ind_sim)  = xn;
end

% Plot quaternion evolution
figure(1),subplot(3,1,1),plot(tsim,xsim(4,:)),grid on, hold on, title('o1')
figure(1),subplot(3,1,2),plot(tsim,xsim(5,:)),grid on, hold on, title('o2')
figure(1),subplot(3,1,3),plot(tsim,xsim(6,:)),grid on, hold on, title('o3')

% Plot chaser angular velocity evolution
figure(2),subplot(3,1,1),plot(tsim,xsim(1,:)),grid on, hold on, title('wcx')
figure(2),subplot(3,1,2),plot(tsim,xsim(2,:)),grid on, hold on, title('wcy')
figure(2),subplot(3,1,3),plot(tsim,xsim(3,:)),grid on, hold on, title('wcz')

close all

%% Discrete time model
n = param.n ;

A = [0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0;
     0, 0, 0, 0, 0, 0;
     1/4, 0, 0, 0, 0, 0;
     0, 1/4, 0, 0, 0, 0;
     0, 0, 1/4, 0, 0, 0];
     
J = param.chaser_inertia;
B = [1/J(1, 1), 0, 0;
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

% system          =   ss(A,B,C,D);
% sys_d           =   c2d(system,Ts,'tustin');
% [Ad,Bd,Cd,Dd]   =   ssdata(sys_d);

[K,S,e] = lqr(A,B,Q,R);

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

xsim        =   zeros(6,Nsim);
xsim(:,1)   =   x0;
tsim        =   0:Ts:(Nsim-1)*Ts;

ulqr        =   zeros(Nsim,3);
uMax        =   0.01;
uMin        =   -0.01;

d           =   [0;0;0] ;

for ind_sim=2:Nsim
    u                   =   - K * (xsim(:,ind_sim-1) - x_ref(:,ind_sim-1));
    u(u>uMax)           =   uMax;
    u(u<uMin)           =   uMin;
%     d                   =   [normrnd(0, 1e-3); normrnd(0, 1e-3); normrnd(0, 1e-3)] ;
    xn                  =   linearAdcsModel(xsim(:,ind_sim-1),u,d,param);
    xsim(:,ind_sim)     =   xn;
    ulqr(ind_sim-1,:)   =   u;
end

% Angle vs reference
figure(4),subplot(3,1,1),plot(tsim,xsim(4,:)),grid on, hold on, title('alfax')
figure(4),subplot(3,1,2),plot(tsim,xsim(5,:)),grid on, hold on, title('alfay')
figure(4),subplot(3,1,3),plot(tsim,xsim(6,:)),grid on, hold on, title('alfaz')

figure(4),subplot(3,1,1),plot(tsim,x_ref(4,1:end-N_MPC),'r--'),grid on, hold on
figure(4),subplot(3,1,2),plot(tsim,x_ref(5,1:end-N_MPC),'r--'),grid on, hold on
figure(4),subplot(3,1,3),plot(tsim,x_ref(6,1:end-N_MPC),'r--'),grid on, hold on

% Angular Velocity vs reference
figure(5),subplot(3,1,1),plot(tsim,xsim(1,:)),grid on, hold on, title('wx')
figure(5),subplot(3,1,2),plot(tsim,xsim(2,:)),grid on, hold on, title('wy')
figure(5),subplot(3,1,3),plot(tsim,xsim(3,:)),grid on, hold on, title('wz')

figure(5),subplot(3,1,1),plot(tsim,x_ref(1,1:end-N_MPC),'r--'),grid on, hold on
figure(5),subplot(3,1,2),plot(tsim,x_ref(2,1:end-N_MPC),'r--'),grid on, hold on
figure(5),subplot(3,1,3),plot(tsim,x_ref(3,1:end-N_MPC),'r--'),grid on, hold on

% close all

figure(6),subplot(2,1,1),plot(tsim,ulqr(:,1)),grid on, hold on, title('Input acceleration [km/s^2] with LQR vs MPC')
figure(6),subplot(2,1,1),plot(tsim,ulqr(:,2)),grid on, hold on
figure(6),subplot(2,1,1),plot(tsim,ulqr(:,3)),grid on, hold on

%% MPC (Close loop simulation w. receding horizon)
Nsim            =   100;
x0              =   ic ;
d               =   [0;0;0];

xMPC            =   zeros(size(ic,1),Nsim);
uMPC            =   zeros(Nsim,3);
xMPC(:,1)       =   x0;

tsim            =   0:Ts:(Nsim-1)*Ts;
ustar           =   [zeros(N,3);zeros(1,3)];
umax            =   [1e-2 1e-2 1e-2] ;

% Set the options and solve the optimization problem in open loop
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
    waitbar(ind_sim/Nsim,f,sprintf('%12.9f',round(ind_sim)))

    ustar               =   fmincon(@(u)linearadc_cost_fun(xMPC(:,ind_sim-1),u,d,N_MPC,Q,R,x_ref,ind_sim,param),[ustar(2:end,:);ustar(end,:)],...
                            [],[],[],[],[],[],@(u)linearadc_constr_fun(xMPC(:,ind_sim-1),u,d,N,umax,param),options);
    u                   =   ustar(1,:)';
%     d                   =   [normrnd(0, 1e-2); normrnd(0, 1e-2); normrnd(0, 1e-2)] ;
    xn                  =   linearAdcsModel(xMPC(:,ind_sim-1),u,d,param);
    xMPC(:,ind_sim)     =   xn;
    uMPC(ind_sim-1,:)   =   u;
end

% Cancel the waitbar when operations finished
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F)

%% Comparison between LQR and MPC trajectories
figure(10),subplot(3,1,1),plot(tsim,xMPC(1,:),'b'),grid on, hold on, title('wx [rad/s]')
figure(10),subplot(3,1,2),plot(tsim,xMPC(2,:),'b'),grid on, hold on, title('wy [rad/s]')
figure(10),subplot(3,1,3),plot(tsim,xMPC(3,:),'b'),grid on, hold on, title('wz [rad/s]')

figure(10),subplot(3,1,1),plot(tsim,xsim(1,:),'k'),grid on, hold on, title('wx [rad/s]')
figure(10),subplot(3,1,2),plot(tsim,xsim(2,:),'k'),grid on, hold on, title('wy [rad/s]')
figure(10),subplot(3,1,3),plot(tsim,xsim(3,:),'k'),grid on, hold on, title('wz [rad/s]')

figure(10),subplot(3,1,1),plot(tsim,x_ref(1,1:end-N_MPC),'r--'),grid on, hold on, title('wx [rad/s]')
figure(10),subplot(3,1,2),plot(tsim,x_ref(2,1:end-N_MPC),'r--'),grid on, hold on, title('wy [rad/s]')
figure(10),subplot(3,1,3),plot(tsim,x_ref(3,1:end-N_MPC),'r--'),grid on, hold on, title('wz [rad/s]')

legend('MPC','LQR','REF')

figure(11),subplot(3,1,1),plot(tsim,xMPC(4,:),'b'),grid on, hold on, title('o1')
figure(11),subplot(3,1,2),plot(tsim,xMPC(5,:),'b'),grid on, hold on, title('o2')
figure(11),subplot(3,1,3),plot(tsim,xMPC(6,:),'b'),grid on, hold on, title('o3')

figure(11),subplot(3,1,1),plot(tsim,xsim(4,:),'k'),grid on, hold on, title('o1')
figure(11),subplot(3,1,2),plot(tsim,xsim(5,:),'k'),grid on, hold on, title('o2')
figure(11),subplot(3,1,3),plot(tsim,xsim(6,:),'k'),grid on, hold on, title('o3')

figure(11),subplot(3,1,1),plot(tsim,x_ref(4,1:end-N_MPC),'r--'),grid on, hold on, title('o1')
figure(11),subplot(3,1,2),plot(tsim,x_ref(5,1:end-N_MPC),'r--'),grid on, hold on, title('o2')
figure(11),subplot(3,1,3),plot(tsim,x_ref(6,1:end-N_MPC),'r--'),grid on, hold on, title('o3')


legend('MPC','LQR','REF')

%% Comparison between FHOCP and MPC control effort
figure(6),subplot(2,1,2),plot(tsim,uMPC(:,1)),grid on, hold on, title('Input acceleration [km/s^2] with LQR vs MPC')
figure(6),subplot(2,1,2),plot(tsim,uMPC(:,2)),grid on, hold on
figure(6),subplot(2,1,2),plot(tsim,uMPC(:,3)),grid on, hold on

legend('LQR', 'MPC')

%% Plot error along the 3-axis
error_x = xMPC(1,:) - x_ref(1,1:end-N_MPC);
error_y = xMPC(2,:) - x_ref(2,1:end-N_MPC);
error_z = xMPC(3,:) - x_ref(3,1:end-N_MPC);

figure,subplot(3,1,1),plot(tsim,error_x,'b','LineWidth', 1.5),grid on
subplot(3,1,2),plot(tsim,error_y,'r','LineWidth', 1.5),grid on
subplot(3,1,3),plot(tsim,error_z,'g','LineWidth', 1.5),grid on

error = [error_x', error_y', error_z'];
norm_err = zeros(100,1);
for i = 1:size(error,1)
    norm_err(i) = norm(error(i,:));
end

figure,plot(tsim,norm_err,'r','LineWidth', 1.5),grid on,title('Norm of the position error')
