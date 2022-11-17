%% Main file for the MPC control of the chaser attitude dynamics.

clc
clear
close all

%% Simulation parameters
param.mu = 398600.0 ;
param.a = 7000 ;
param.chaser_inertia = [0.872, 0, 0;
                        0, 0.115, 0;
                        0, 0, 0.797] ;
param.target_inertia = [0.5208, 0, 0;
                        0, 0.5208, 0;
                        0, 0, 0.6667] ;
param.wmax = 0.1 ; % [rad/s]
                    
% Sampling time
Ts = 1 ; % [s]

% Initial conditions (chaser attitude + target attitude)
ic = [0; 0; 0; 1; 0; 0; 0; deg2rad(50); deg2rad(10); deg2rad(-30); 0.01; 0.002; -0.005] ;

%% Open loop simulation
Nsim        =   30 ;

x0          =   ic ;

u           =   [1e-3; 0; 0] ;
d           =   [0.0; -0.0; 0.0] ;

xsim        =   zeros(size(ic, 1), Nsim) ;
xsim(:,1)   =   x0 ;

tsim        =   0:Ts:(Nsim-1)*Ts ;

for ind_sim=2:Nsim
    [tout, xout]     = ode45(@(t,x)adcmodel(t,x,u,d,param), [0 Ts], xsim(:,ind_sim-1)) ;
    xsim(:,ind_sim)  = xout(end,:)' ;
end

% Plot quaternion evolution
figure(1),subplot(5,1,1),plot(tsim,xsim(1,:)),grid on, hold on, title('q1')
figure(1),subplot(5,1,2),plot(tsim,xsim(2,:)),grid on, hold on, title('q2')
figure(1),subplot(5,1,3),plot(tsim,xsim(3,:)),grid on, hold on, title('q3')
figure(1),subplot(5,1,4),plot(tsim,xsim(4,:)),grid on, hold on, title('q4')
figure(1),subplot(5,1,5),plot(tsim,u*ones(1,Nsim)),grid on, hold on, title('Input Control Torque [Nm]')

% Plot chaser angular velocity evolution
figure(2),subplot(3,1,1),plot(tsim,xsim(5,:)),grid on, hold on, title('wcx')
figure(2),subplot(3,1,2),plot(tsim,xsim(6,:)),grid on, hold on, title('wcy')
figure(2),subplot(3,1,3),plot(tsim,xsim(7,:)),grid on, hold on, title('wcz')

% Plot target euler angles evolution
figure(3),subplot(3,1,1),plot(tsim,xsim(8,:)), grid on, hold on, title('alfax')
figure(3),subplot(3,1,2),plot(tsim,xsim(9,:)), grid on, hold on, title('alfay')
figure(3),subplot(3,1,3),plot(tsim,xsim(10,:)),grid on, hold on, title('alfaz')

close all

%% Finite Horizon Optimal Control Problem
x0          =   ic ;
N           =   40 ;
Q           =   1 ;
R           =   1e-5 ;
umax        =   [1e-2 1e-2 1e-2] ;
x_r         =   ic(8:end) ;
% Transform euler angles to quaternion
x_ref_quat  =   DCM_quat(Angles321_DCM(x_r(1:3))) ;

x_ref       =   zeros(size(x_r,1)+1,N) ;
x_ref(:,1)  =   [x_ref_quat; x_r(4:6)] ;

options     =   optimset('Display','Iter','MaxFunEvals',1e4,'Algorithm','active-set');  

ustar   = fmincon(@(u)adc_cost_fun(x0,u,d,N,Q,R,x_ref(:,1),Ts,param),zeros(N,3),...
                    [],[],[],[],[],[],@(u)adc_constr_fun(x0,u,d,N,umax,Ts,param),options);

%figure,plot(ustar),grid on

[Constr,Constr_eq,xpred]=adc_constr_fun(x0,ustar,d,N,umax,Ts,param);

%% Open loop simulation - FHOCP solution
Nsim            =   N;
%d               =   [0; 0; 0];

xsim            =   zeros(size(x0,1),Nsim);
xsim(:,1)       =   x0;
tsim            =   0:Ts:(Nsim-1)*Ts;

for ind_sim=2:Nsim
    u                = ustar(ind_sim-1,:)';
    [tout,xout]      = ode45(@(t,x)adcmodel(t,x,u,d,param),[0 Ts],xsim(:,ind_sim-1));
    xsim(:,ind_sim)  = xout(end,:)';
    x_r              = xsim(8:end, ind_sim);
    x_ref_quat       = DCM_quat(Angles321_DCM(x_r(1:3))) ;
    x_ref(:,ind_sim) = [x_ref_quat; x_r(4:6)] ;
end

% Plot quaternion evolution with FHOCP
figure(1),subplot(5,1,1),plot(tsim,xsim(1,:)),grid on, hold on, title('q1')
figure(1),subplot(5,1,2),plot(tsim,xsim(2,:)),grid on, hold on, title('q2')
figure(1),subplot(5,1,3),plot(tsim,xsim(3,:)),grid on, hold on, title('q3')
figure(1),subplot(5,1,4),plot(tsim,xsim(4,:)),grid on, hold on, title('q4')
figure(1),subplot(5,1,5),plot(tsim,ustar*ones(3,Nsim)),grid on, hold on, title('Input Control Torque [Nm]')

% Plot chaser angular velocity evolution with FHOCP
figure(2),subplot(3,1,1),plot(tsim,xsim(5,:)),grid on, hold on, title('wcx')
figure(2),subplot(3,1,2),plot(tsim,xsim(6,:)),grid on, hold on, title('wcy')
figure(2),subplot(3,1,3),plot(tsim,xsim(7,:)),grid on, hold on, title('wcz')

% Plot target euler angles evolution with FHOCP
% figure(3),subplot(3,1,1),plot(tsim,xsim(8,:)), grid on, hold on, title('alfax')
% figure(3),subplot(3,1,2),plot(tsim,xsim(9,:)), grid on, hold on, title('alfay')
% figure(3),subplot(3,1,3),plot(tsim,xsim(10,:)),grid on, hold on, title('alfaz')

figure(4),subplot(2,1,1),plot(tsim,ustar*ones(3,Nsim)),grid on, hold on, title('Input Control Torque with FHOCP [Nm]')

% Plot quaternion evolution against reference
figure(5),subplot(4,1,1),plot(tsim,xsim(1,:)),grid on, hold on, title('q1 vs ref')
figure(5),subplot(4,1,2),plot(tsim,xsim(2,:)),grid on, hold on, title('q2 vs ref')
figure(5),subplot(4,1,3),plot(tsim,xsim(3,:)),grid on, hold on, title('q3 vs ref')
figure(5),subplot(4,1,4),plot(tsim,xsim(4,:)),grid on, hold on, title('q4 vs ref')
figure(5),subplot(4,1,1),plot(tsim,x_ref(1,:),'r--'),grid on, hold on
figure(5),subplot(4,1,2),plot(tsim,x_ref(2,:),'r--'),grid on, hold on
figure(5),subplot(4,1,3),plot(tsim,x_ref(3,:),'r--'),grid on, hold on
figure(5),subplot(4,1,4),plot(tsim,x_ref(4,:),'r--'),grid on, hold on

%% MPC (Close loop simulation w. receding horizon)
Nsim        =   N;
N_MPC       =   30;
x0          =   ic ;
%d           =   [0; 0; 0];

xMPC        =   zeros(size(x0,1),Nsim);
uMPC        =   zeros(Nsim,3);
xMPC(:,1)   =   x0;
tsim        =   0:Ts:(Nsim-1)*Ts;
ustar       =   [zeros(N,3);zeros(1,3)];

x_ref       =   zeros(size(x_r,1)+1,N) ;
x_ref(:,1)  =   [x_ref_quat; x_r(4:6)] ;

for ind_sim=2:Nsim
    ustar               =   fmincon(@(u)adc_cost_fun(xMPC(:,ind_sim-1),u,d,N_MPC,Q,R,x_ref(:,ind_sim-1),Ts,param),[ustar(2:end,:);ustar(end,:)],...
                                    [],[],[],[],[],[],@(u)adc_constr_fun(xMPC(:,ind_sim-1),u,d,N_MPC,umax,Ts,param),options);
    u                   =   ustar(1,:)';
    [tout,xout]         =   ode45(@(t,x)adcmodel(t,x,u,d,param),[0 Ts],xMPC(:,ind_sim-1));
    xMPC(:,ind_sim)     =   xout(end,:)';
    x_r                 =   xsim(8:end, ind_sim);
    x_ref_quat          =   DCM_quat(Angles321_DCM(x_r(1:3))) ;
    x_ref(:,ind_sim)    =   [x_ref_quat; x_r(4:6)] ;
    uMPC(ind_sim-1,:)   =   u;
end
%%
% Plot quaternion evolution with MPC
figure(1),subplot(5,1,1),plot(tsim,xMPC(1,:)),grid on, hold on, title('q1')
figure(1),subplot(5,1,2),plot(tsim,xMPC(2,:)),grid on, hold on, title('q2')
figure(1),subplot(5,1,3),plot(tsim,xMPC(3,:)),grid on, hold on, title('q3')
figure(1),subplot(5,1,4),plot(tsim,xMPC(4,:)),grid on, hold on, title('q4')
figure(1),subplot(5,1,5),plot(tsim,uMPC*ones(3,Nsim)),grid on, hold on, title('Input Control Torque [Nm]')

% Plot chaser angular velocity evolution with MPC
figure(2),subplot(3,1,1),plot(tsim,xMPC(5,:)),grid on, hold on, title('wcx')
figure(2),subplot(3,1,2),plot(tsim,xMPC(6,:)),grid on, hold on, title('wcy')
figure(2),subplot(3,1,3),plot(tsim,xMPC(7,:)),grid on, hold on, title('wcz')

% Plot target euler angles evolution with MPC
% figure(3),subplot(3,1,1),plot(tsim,xMPC(8,:)), grid on, hold on, title('alfax')
% figure(3),subplot(3,1,2),plot(tsim,xMPC(9,:)), grid on, hold on, title('alfay')
% figure(3),subplot(3,1,3),plot(tsim,xMPC(10,:)),grid on, hold on, title('alfaz')

% Plot quaternion evolution against reference
figure(6),subplot(4,1,1),plot(tsim,xMPC(1,:)),grid on, hold on, title('q1 vs ref')
figure(6),subplot(4,1,2),plot(tsim,xMPC(2,:)),grid on, hold on, title('q2 vs ref')
figure(6),subplot(4,1,3),plot(tsim,xMPC(3,:)),grid on, hold on, title('q3 vs ref')
figure(6),subplot(4,1,4),plot(tsim,xMPC(4,:)),grid on, hold on, title('q4 vs ref')
figure(6),subplot(4,1,1),plot(tsim,x_ref(1,:),'r--'),grid on, hold on
figure(6),subplot(4,1,2),plot(tsim,x_ref(2,:),'r--'),grid on, hold on
figure(6),subplot(4,1,3),plot(tsim,x_ref(3,:),'r--'),grid on, hold on
figure(6),subplot(4,1,4),plot(tsim,x_ref(4,:),'r--'),grid on, hold on

%% Control comparison between FHOCP and MPC

figure(4),subplot(2,1,2),plot(tsim,uMPC*ones(3,Nsim)),grid on, hold on, title('Input Control Torque with MPC [Nm]')
