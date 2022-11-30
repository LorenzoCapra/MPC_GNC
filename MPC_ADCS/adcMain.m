%% Main file for the MPC control of the chaser attitude dynamics.

clc
clear
close all

%% Simulation parameters
param.mu = 398600.0 ;
param.Re = 6378.0 ;                     % Earth radius
param.h = 500.0 ;                       % Orbit altitude
param.R = param.Re + param.h ;          % Orbit radius
param.n = sqrt(param.mu/(param.R^3)) ;  % Mean orbital motion
param.chaser_inertia = [0.872, 0, 0;
                        0, 0.115, 0;
                        0, 0, 0.797] ;
param.target_inertia = [0.5208, 0, 0;
                        0, 0.5208, 0;
                        0, 0, 0.6667] ;
param.wmax = 0.1 ; % [rad/s]
                    
% Sampling time
Ts = 1 ; % [s]

% Initial conditions:
alfa        =   [0; 0; 0];
w_target    =   [0.0; 0.0; 0.05];
target_ic   =   [alfa; w_target];

q0          =   [0; 0; 0; 1];
w_chaser    =   [0; 0; 0];

ic          =   [q0; w_chaser] ;

%% Open loop simulation
Nsim        =   30 ;

x0          =   ic ;

u           =   [normrnd(0, 1e-3); normrnd(0, 1e-3); normrnd(0, 1e-3)] ;
d           =   [normrnd(0, 1e-4); normrnd(0, 1e-4); normrnd(0, 1e-4)] ;

xsim        =   zeros(size(ic, 1), Nsim) ;
xsim(:,1)   =   x0 ;

tsim        =   0:Ts:(Nsim-1)*Ts ;

for ind_sim=2:Nsim
    [tout, xout]     = ode45(@(t,x)adcmodel(t,x,u,d,param), [0 Ts], xsim(:,ind_sim-1)) ;
    xsim(:,ind_sim)  = xout(end,:)' ;
end

% Plot quaternion evolution
figure(1),subplot(4,1,1),plot(tsim,xsim(1,:)),grid on, hold on, title('q1')
figure(1),subplot(4,1,2),plot(tsim,xsim(2,:)),grid on, hold on, title('q2')
figure(1),subplot(4,1,3),plot(tsim,xsim(3,:)),grid on, hold on, title('q3')
figure(1),subplot(4,1,4),plot(tsim,xsim(4,:)),grid on, hold on, title('q4')

% Plot chaser angular velocity evolution
figure(2),subplot(3,1,1),plot(tsim,xsim(5,:)),grid on, hold on, title('wcx')
figure(2),subplot(3,1,2),plot(tsim,xsim(6,:)),grid on, hold on, title('wcy')
figure(2),subplot(3,1,3),plot(tsim,xsim(7,:)),grid on, hold on, title('wcz')


%% Reference attitude (qauternion + angular velocity)
N           =   40 ;
N_MPC       =   10;

Tt          =   0:Ts:Ts*(N-1+N_MPC);
[~,xtout]   =   ode45(@(t,x)targetode(t,x,param),Tt,target_ic);
x_r         =   xtout';

% Transform euler angles to quaternion
x_ref         =   zeros(size(ic,1),N+N_MPC);
for i = 1:N+N_MPC
    x_ref_quat  =   DCM_quat(Angles321_DCM(x_r(1:3,i))) ;
    x_ref(:,i)  =   [x_ref_quat; x_r(4:6,i)];
end

% Plot target quaternion evolution
figure(3),subplot(4,1,1),plot(Tt,x_ref(1,:)), grid on, hold on, title('q1 ref')
figure(3),subplot(4,1,2),plot(Tt,x_ref(2,:)), grid on, hold on, title('q2 ref')
figure(3),subplot(4,1,3),plot(Tt,x_ref(3,:)), grid on, hold on, title('q3 ref')
figure(3),subplot(4,1,4),plot(Tt,x_ref(4,:)), grid on, hold on, title('q4 ref')

% Plot target angular velocity evolution
figure(4),subplot(3,1,1),plot(Tt,x_ref(5,:)), grid on, hold on, title('wx ref')
figure(4),subplot(3,1,2),plot(Tt,x_ref(6,:)), grid on, hold on, title('wy ref')
figure(4),subplot(3,1,3),plot(Tt,x_ref(7,:)), grid on, hold on, title('wz ref')

%% Finite Horizon Optimal Control Problem
close all

x0          =   ic ;
Q           =   1 ;
R           =   1e-5 ;
umax        =   [1e-2 1e-2 1e-2] ;

options     =   optimset('Display','Iter','MaxFunEvals',1e4,'Algorithm','active-set');  

ustar   = fmincon(@(u)adc_cost_fun(x0,u,d,N,Q,R,x_ref,Ts,param,1),zeros(N,3),...
                    [],[],[],[],[],[],@(u)adc_constr_fun(x0,u,d,N,umax,Ts,param),options);

%% Open loop simulation - FHOCP solution
Nsim            =   N;
%d              =   [0; 0; 0];

xsim            =   zeros(size(x0,1),Nsim);
xsim(:,1)       =   x0;
tsim            =   0:Ts:(Nsim-1)*Ts;

for ind_sim=2:Nsim
    u                = ustar(ind_sim-1,:)';
    [tout,xout]      = ode45(@(t,x)adcmodel(t,x,u,d,param),[0 Ts],xsim(:,ind_sim-1));
    xsim(:,ind_sim)  = xout(end,:)';
end

% Plot quaternion evolution with FHOCP
figure(1),subplot(5,1,1),plot(tsim,xsim(1,:)),grid on, hold on, title('q1')
figure(1),subplot(5,1,2),plot(tsim,xsim(2,:)),grid on, hold on, title('q2')
figure(1),subplot(5,1,3),plot(tsim,xsim(3,:)),grid on, hold on, title('q3')
figure(1),subplot(5,1,4),plot(tsim,xsim(4,:)),grid on, hold on, title('q4')

% Plot chaser angular velocity evolution with FHOCP
figure(2),subplot(3,1,1),plot(tsim,xsim(5,:)),grid on, hold on, title('wcx')
figure(2),subplot(3,1,2),plot(tsim,xsim(6,:)),grid on, hold on, title('wcy')
figure(2),subplot(3,1,3),plot(tsim,xsim(7,:)),grid on, hold on, title('wcz')

% Plot target euler angles evolution with FHOCP
% figure(3),subplot(3,1,1),plot(tsim,xsim(8,:)), grid on, hold on, title('alfax')
% figure(3),subplot(3,1,2),plot(tsim,xsim(9,:)), grid on, hold on, title('alfay')
% figure(3),subplot(3,1,3),plot(tsim,xsim(10,:)),grid on, hold on, title('alfaz')

% Plot quaternion evolution against reference
figure(5),subplot(4,1,1),plot(tsim,xsim(1,:),'k'),grid on, hold on, title('q1 vs ref')
figure(5),subplot(4,1,2),plot(tsim,xsim(2,:),'k'),grid on, hold on, title('q2 vs ref')
figure(5),subplot(4,1,3),plot(tsim,xsim(3,:),'k'),grid on, hold on, title('q3 vs ref')
figure(5),subplot(4,1,4),plot(tsim,xsim(4,:),'k'),grid on, hold on, title('q4 vs ref')
figure(5),subplot(4,1,1),plot(tsim,x_ref(1,1:end-N_MPC),'r--'),grid on, hold on
figure(5),subplot(4,1,2),plot(tsim,x_ref(2,1:end-N_MPC),'r--'),grid on, hold on
figure(5),subplot(4,1,3),plot(tsim,x_ref(3,1:end-N_MPC),'r--'),grid on, hold on
figure(5),subplot(4,1,4),plot(tsim,x_ref(4,1:end-N_MPC),'r--'),grid on, hold on

% Plot angular velocity evolution against reference
figure(7),subplot(3,1,1),plot(tsim,xsim(5,:),'k'),grid on, hold on, title('wx vs ref')
figure(7),subplot(3,1,2),plot(tsim,xsim(6,:),'k'),grid on, hold on, title('wy vs ref')
figure(7),subplot(3,1,3),plot(tsim,xsim(7,:),'k'),grid on, hold on, title('wz vs ref')
figure(7),subplot(3,1,1),plot(tsim,x_ref(5,1:end-N_MPC),'r--'),grid on, hold on
figure(7),subplot(3,1,2),plot(tsim,x_ref(6,1:end-N_MPC),'r--'),grid on, hold on
figure(7),subplot(3,1,3),plot(tsim,x_ref(7,1:end-N_MPC),'r--'),grid on, hold on

% Plot control effort
figure(666),subplot(3,1,1),plot(tsim,ustar(:,1)),grid on, hold on, title('Input acceleration [m/s^2] with MPC vs FHOCP')
figure(666),subplot(3,1,2),plot(tsim,ustar(:,2)),grid on, hold on
figure(666),subplot(3,1,3),plot(tsim,ustar(:,3)),grid on, hold on

%% MPC (Close loop simulation w. receding horizon)
Nsim        =   N;
N_MPC       =   10;
x0          =   ic ;
%d           =   [0; 0; 0];

xMPC        =   zeros(size(x0,1),Nsim);
uMPC        =   zeros(Nsim,3);
xMPC(:,1)   =   x0;
tsim        =   0:Ts:(Nsim-1)*Ts;
ustar       =   [zeros(N,3);zeros(1,3)];

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

    ustar               =   fmincon(@(u)adc_cost_fun(xMPC(:,ind_sim-1),u,d,N_MPC,Q,R,x_ref,Ts,param,ind_sim),[ustar(2:end,:);ustar(end,:)],...
                                    [],[],[],[],[],[],@(u)adc_constr_fun(xMPC(:,ind_sim-1),u,d,N_MPC,umax,Ts,param),options);
    u                   =   ustar(1,:)';
    [tout,xout]         =   ode45(@(t,x)adcmodel(t,x,u,d,param),[0 Ts],xMPC(:,ind_sim-1));
    xMPC(:,ind_sim)     =   xout(end,:)';
    uMPC(ind_sim-1,:)   =   u;
end

delete(f)

%% Plot quaternion evolution with MPC
figure(1),subplot(5,1,1),plot(tsim,xMPC(1,:)),grid on, hold on, title('q1')
figure(1),subplot(5,1,2),plot(tsim,xMPC(2,:)),grid on, hold on, title('q2')
figure(1),subplot(5,1,3),plot(tsim,xMPC(3,:)),grid on, hold on, title('q3')
figure(1),subplot(5,1,4),plot(tsim,xMPC(4,:)),grid on, hold on, title('q4')

% Plot chaser angular velocity evolution with MPC
figure(2),subplot(3,1,1),plot(tsim,xMPC(5,:)),grid on, hold on, title('wcx')
figure(2),subplot(3,1,2),plot(tsim,xMPC(6,:)),grid on, hold on, title('wcy')
figure(2),subplot(3,1,3),plot(tsim,xMPC(7,:)),grid on, hold on, title('wcz')

% Plot target euler angles evolution with MPC
% figure(3),subplot(3,1,1),plot(tsim,xMPC(8,:)), grid on, hold on, title('alfax')
% figure(3),subplot(3,1,2),plot(tsim,xMPC(9,:)), grid on, hold on, title('alfay')
% figure(3),subplot(3,1,3),plot(tsim,xMPC(10,:)),grid on, hold on, title('alfaz')

% Plot quaternion evolution against reference
figure(5),subplot(4,1,1),plot(tsim,xMPC(1,:),'b'),grid on, hold on, title('q1 vs ref')
figure(5),subplot(4,1,2),plot(tsim,xMPC(2,:),'b'),grid on, hold on, title('q2 vs ref')
figure(5),subplot(4,1,3),plot(tsim,xMPC(3,:),'b'),grid on, hold on, title('q3 vs ref')
figure(5),subplot(4,1,4),plot(tsim,xMPC(4,:),'b'),grid on, hold on, title('q4 vs ref')

% Plot angular velocity evolution against reference
figure(7),subplot(3,1,1),plot(tsim,xMPC(5,:),'b'),grid on, hold on,
figure(7),subplot(3,1,2),plot(tsim,xMPC(6,:),'b'),grid on, hold on,
figure(7),subplot(3,1,3),plot(tsim,xMPC(7,:),'b'),grid on, hold on,

%% Comparison between FHOCP and MPC control effort
figure(666),subplot(3,1,1),plot(tsim,uMPC(:,1)),grid on, hold on, title('Input acceleration [m/s^2] with MPC vs FHOCP')
figure(666),subplot(3,1,2),plot(tsim,uMPC(:,2)),grid on, hold on
figure(666),subplot(3,1,3),plot(tsim,uMPC(:,3)),grid on, hold on

%% Plot error along the 3-axis
error_q1        = xMPC(1,:) - x_ref(1,1:end-N_MPC);
error_q2        = xMPC(2,:) - x_ref(2,1:end-N_MPC);
error_q3        = xMPC(3,:) - x_ref(3,1:end-N_MPC);
error_q4        = xMPC(4,:) - x_ref(4,1:end-N_MPC);
error_wx        = xMPC(5,:) - x_ref(5,1:end-N_MPC);
error_wy        = xMPC(6,:) - x_ref(6,1:end-N_MPC);
error_wz        = xMPC(7,:) - x_ref(7,1:end-N_MPC);

e_q1            = xsim(1,:) - x_ref(1,1:end-N_MPC);
e_q2            = xsim(2,:) - x_ref(2,1:end-N_MPC);
e_q3            = xsim(3,:) - x_ref(3,1:end-N_MPC);
e_q4            = xsim(4,:) - x_ref(4,1:end-N_MPC);
e_wx            = xsim(5,:) - x_ref(5,1:end-N_MPC);
e_wy            = xsim(6,:) - x_ref(6,1:end-N_MPC);
e_wz            = xsim(7,:) - x_ref(7,1:end-N_MPC);

figure(8),subplot(4,1,1),plot(tsim,error_q1,'b'),grid on,hold on,title('Quaternion error (MPC vs FHOCP)')
figure(8),subplot(4,1,2),plot(tsim,error_q2,'r'),grid on,hold on
figure(8),subplot(4,1,3),plot(tsim,error_q3,'g'),grid on,hold on
figure(8),subplot(4,1,4),plot(tsim,error_q4,'m'),grid on,hold on

figure(8),subplot(4,1,1),plot(tsim,e_q1,'b--'),grid on,hold on
figure(8),subplot(4,1,2),plot(tsim,e_q2,'r--'),grid on,hold on
figure(8),subplot(4,1,3),plot(tsim,e_q3,'g--'),grid on,hold on
figure(8),subplot(4,1,4),plot(tsim,e_q4,'m--'),grid on,hold on

error_angleMPC   = [error_q1', error_q2', error_q3', error_q4'];
error_angleFHOCP = [e_q1', e_q2', e_q3', e_q4'];
error_wMPC   = [error_wx', error_wy', error_wz'];
error_wFHOCP = [e_wx', e_wy', e_wz'];
norm_err = zeros(N,4);
for i = 1:N
    norm_err(i,1) = norm(error_angleMPC(i,:));
    norm_err(i,2) = norm(error_angleFHOCP(i,:));
    norm_err(i,3) = norm(error_wMPC(i,:));
    norm_err(i,4) = norm(error_wFHOCP(i,:));
end

figure(16),subplot(2,1,1),plot(tsim,norm_err(:,1),'r'),grid on,hold on,title('Norm of the quaternion error (MPC vs FHOCP)')
figure(16),subplot(2,1,1),plot(tsim,norm_err(:,2),'b'),grid on,hold on
figure(16),subplot(2,1,2),plot(tsim,norm_err(:,3),'r'),grid on,hold on,title('Norm of the angular velocity error (MPC vs FHOCP)')
figure(16),subplot(2,1,2),plot(tsim,norm_err(:,4),'b'),grid on,hold on
