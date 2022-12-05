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

param.Jinv = inv(param.chaser_inertia);
                    
% Sampling time
Ts = 1 ; % [s]

% Initial conditions:
q0          =   [0; 0; 0; 1];
w_chaser    =   [0; 0; 0];

ic          =   [q0(1:3); w_chaser] ;

%% Open loop simulation
Nsim        =   50 ;

x0          =   ic ;

u           =   [normrnd(0, 1e-3); normrnd(0, 1e-3); normrnd(0, 1e-3)] ;
d           =   [normrnd(0, 1e-4); normrnd(0, 1e-4); normrnd(0, 1e-4)] ;

xsim        =   zeros(size(ic, 1), Nsim) ;
xsim(:,1)   =   x0 ;
q4          =   zeros(1,Nsim);
q4(:,1)     =   sqrt(1 - (x0(1)^2) - (x0(2)^2) - (x0(3)^2));

tsim        =   0:Ts:(Nsim-1)*Ts ;

for ind_sim=2:Nsim
    [tout, xout]     = ode45(@(t,x)linearAdcsModel(t,x,u,d,param), [0 Ts], xsim(:,ind_sim-1)) ;
    xs               = xout(end,:)' ;
    xsim(:,ind_sim)  = xs;
    q4(:,ind_sim)    = sqrt(1 - (xs(1)^2) - (xs(2)^2) - (xs(3)^2));
end

% Plot quaternion evolution
figure(1),subplot(4,1,1),plot(tsim,xsim(1,:)),grid on, hold on, title('q1')
figure(1),subplot(4,1,2),plot(tsim,xsim(2,:)),grid on, hold on, title('q2')
figure(1),subplot(4,1,3),plot(tsim,xsim(3,:)),grid on, hold on, title('q3')
figure(1),subplot(4,1,4),plot(tsim,q4(1,:)),grid on, hold on, title('q4')

% Plot chaser angular velocity evolution
figure(2),subplot(3,1,1),plot(tsim,xsim(4,:)),grid on, hold on, title('wcx')
figure(2),subplot(3,1,2),plot(tsim,xsim(5,:)),grid on, hold on, title('wcy')
figure(2),subplot(3,1,3),plot(tsim,xsim(6,:)),grid on, hold on, title('wcz')

close all

%% Discretize the state space model
A = [0.5*eye(3,3) , zeros(3,3);
     zeros(3,3)   , zeros(3,3)] ;
 
Jinv = param.Jinv;
B = [zeros(3,3); Jinv];

C = eye(6);

D = 0 ;

system          =   ss(A,B,C,D);
sys_d           =   c2d(system,Ts,'tustin');
[Ad,Bd,Cd,Dd]   =   ssdata(sys_d);

%% Reference attitude (qauternion + angular velocity)
N           =   Nsim+10 ;
N_MPC       =   10;
Tt          =   0:Ts:Ts*(N-1+N_MPC);

qt          =   [0; 0; 0; 1];
w_target    =   [0.0; 0.0; 0.05];
target_ic   =   [qt; w_target];

x_ref         =     zeros(size(ic,1),N+N_MPC);
[~,xtout]     =     ode45(@(t,x)targetode(t,x,param),Tt,target_ic);
x_ref(1:7,:)  =     xtout';

% Plot target quaternion evolution
figure(3),subplot(4,1,1),plot(Tt,x_ref(1,:)), grid on, hold on, title('q1 ref')
figure(3),subplot(4,1,2),plot(Tt,x_ref(2,:)), grid on, hold on, title('q2 ref')
figure(3),subplot(4,1,3),plot(Tt,x_ref(3,:)), grid on, hold on, title('q3 ref')
figure(3),subplot(4,1,4),plot(Tt,x_ref(4,:)), grid on, hold on, title('q4 ref')

% Plot target angular velocity evolution
figure(4),subplot(3,1,1),plot(Tt,x_ref(5,:)), grid on, hold on, title('wx ref')
figure(4),subplot(3,1,2),plot(Tt,x_ref(6,:)), grid on, hold on, title('wy ref')
figure(4),subplot(3,1,3),plot(Tt,x_ref(7,:)), grid on, hold on, title('wz ref')

% Remove q4 from reference array
x_ref_c = x_ref;
x_ref(4,:) = [];

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x_ref = x_ref(1:6, 40);    %%
% x_ref(end) = 0;            %%
% x_ref_c = x_ref_c(4, 40);  %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

%% Finite Horizon Optimal Control Problem
x0          =   ic ;
N           =   60 ;
Q           =   1;
R           =   1e-5 ;
umax        =   [1e-2 1e-2 1e-2] ;
d           =   [0;0;0];

% Set the options and solve the optimization problem in open loop
options     =   optimset('Display','Iter','MaxFunEvals',1e4,'Algorithm','active-set'); 

ustar   = fmincon(@(u)linearadc_cost_fun(x0,u,d,N,Q,R,x_ref,1,Ad,Bd),zeros(N,3),...
                [],[],[],[],[],[],@(u)linearadc_constr_fun(x0,u,N,umax,param,Ad,Bd),options);
            
%% Open loop simulation - FHOCP solution
Nsim            =   N;
d               =   [0; 0; 0];

xsim            =   zeros(size(x0,1),Nsim);
xsim(:,1)       =   x0;
q4_sim          =   zeros(1,Nsim);
q4_sim(:,1)     =   sqrt(1 - (x0(1)^2) - (x0(2)^2) - (x0(3)^2));
tsim            =   0:Ts:(Nsim-1)*Ts;

for ind_sim=2:Nsim
    u                   = ustar(ind_sim-1,:)';
    [tout,xout]         = ode45(@(t,x)linearAdcsModel(t,x,u,d,param),[0 Ts],xsim(:,ind_sim-1));
    xs                  = xout(end,:)';
    xsim(:,ind_sim)     = xs;
    q4_sim(:,ind_sim)   = sqrt(1 - (xs(1)^2) - (xs(2)^2) - (xs(3)^2));
end

% Plot quaternion evolution with FHOCP
figure(1),subplot(4,1,1),plot(tsim,xsim(1,:)),grid on, hold on, title('q1')
figure(1),subplot(4,1,2),plot(tsim,xsim(2,:)),grid on, hold on, title('q2')
figure(1),subplot(4,1,3),plot(tsim,xsim(3,:)),grid on, hold on, title('q3')
figure(1),subplot(4,1,4),plot(tsim,q4_sim(1,:)),grid on, hold on, title('q4')

% Plot chaser angular velocity evolution with FHOCP
figure(2),subplot(3,1,1),plot(tsim,xsim(4,:)),grid on, hold on, title('wcx')
figure(2),subplot(3,1,2),plot(tsim,xsim(5,:)),grid on, hold on, title('wcy')
figure(2),subplot(3,1,3),plot(tsim,xsim(6,:)),grid on, hold on, title('wcz')

% Plot target euler angles evolution
% figure(1),subplot(3,1,1),plot(tsim,xsim(8,:)), grid on, hold on, title('alfax')
% figure(1),subplot(3,1,2),plot(tsim,xsim(9,:)), grid on, hold on, title('alfay')
% figure(1),subplot(3,1,3),plot(tsim,xsim(10,:)),grid on, hold on, title('alfaz')

% Plot quaternion evolution against reference
figure(5),subplot(4,1,1),plot(tsim,xsim(1,:),'k'),grid on, hold on, title('q1 vs ref')
figure(5),subplot(4,1,2),plot(tsim,xsim(2,:),'k'),grid on, hold on, title('q2 vs ref')
figure(5),subplot(4,1,3),plot(tsim,xsim(3,:),'k'),grid on, hold on, title('q3 vs ref')
figure(5),subplot(4,1,4),plot(tsim,q4_sim(1,:),'k'),grid on, hold on, title('q4 vs ref')
% figure(5),subplot(4,1,1),yline(x_ref(1),'r--'),grid on, hold on,
% figure(5),subplot(4,1,2),yline(x_ref(2),'r--'),grid on, hold on,
% figure(5),subplot(4,1,3),yline(x_ref(3),'r--'),grid on, hold on,
% figure(5),subplot(4,1,4),yline(x_ref_c,'r--'),grid on, hold on,
figure(5),subplot(4,1,1),plot(tsim,x_ref(1,1:end-N_MPC),'r--'),grid on, hold on
figure(5),subplot(4,1,2),plot(tsim,x_ref(2,1:end-N_MPC),'r--'),grid on, hold on
figure(5),subplot(4,1,3),plot(tsim,x_ref(3,1:end-N_MPC),'r--'),grid on, hold on
figure(5),subplot(4,1,4),plot(tsim,x_ref_c(4,1:end-N_MPC),'r--'),grid on, hold on

% Plot angular velocity evolution against reference
figure(7),subplot(3,1,1),plot(tsim,xsim(4,:),'k'),grid on, hold on, title('wx vs ref')
figure(7),subplot(3,1,2),plot(tsim,xsim(5,:),'k'),grid on, hold on, title('wy vs ref')
figure(7),subplot(3,1,3),plot(tsim,xsim(6,:),'k'),grid on, hold on, title('wz vs ref')
% figure(7),subplot(3,1,1),yline(x_ref(4),'r--'),grid on, hold on,
% figure(7),subplot(3,1,2),yline(x_ref(5),'r--'),grid on, hold on,
% figure(7),subplot(3,1,3),yline(x_ref(6),'r--'),grid on, hold on,
figure(7),subplot(3,1,1),plot(tsim,x_ref(4,1:end-N_MPC),'r--'),grid on, hold on
figure(7),subplot(3,1,2),plot(tsim,x_ref(5,1:end-N_MPC),'r--'),grid on, hold on
figure(7),subplot(3,1,3),plot(tsim,x_ref(6,1:end-N_MPC),'r--'),grid on, hold on

close all

% Plot control effort
figure(666),subplot(3,1,1),plot(tsim,ustar(:,1)),grid on, hold on, title('Input acceleration [m/s^2] with MPC vs FHOCP')
figure(666),subplot(3,1,2),plot(tsim,ustar(:,2)),grid on, hold on
figure(666),subplot(3,1,3),plot(tsim,ustar(:,3)),grid on, hold on

%% MPC (Close loop simulation w. receding horizon)
Nsim            =   N;
x0              =   ic ;
d               =   [0;0;0];

xMPC            =   zeros(size(ic,1),Nsim);
uMPC            =   zeros(Nsim,3);
xMPC(:,1)       =   x0;
q4_mpc          =   zeros(1,Nsim);
q4_mpc(:,1)     =   sqrt(1 - (x0(1)^2) - (x0(2)^2) - (x0(3)^2));
tsim            =   0:Ts:(Nsim-1)*Ts;
ustar           =   [zeros(N,3);zeros(1,3)];

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

    ustar               =   fmincon(@(u)linearadc_cost_fun(xMPC(:,ind_sim-1),u,d,N_MPC,Q,R,x_ref,ind_sim,A,B),[ustar(2:end,:);ustar(end,:)],...
                            [],[],[],[],[],[],@(u)linearadc_constr_fun(xMPC(:,ind_sim-1),u,N,umax,param,A,B),options);
    u                   =   ustar(1,:)';
    [tout,xout]         =   ode45(@(t,x)linearAdcsModel(t,x,u,d,param),[0 Ts],xMPC(:,ind_sim-1));
    xs                  =   xout(end,:)';
    xMPC(:,ind_sim)     =   xs;
    q4_mpc(:,ind_sim)   =   sqrt(1 - (xs(1)^2) - (xs(2)^2) - (xs(3)^2));
    uMPC(ind_sim-1,:)   =   u;
end

% Cancel the waitbar when operations finished
F = findall(0,'type','figure','tag','TMWWaitbar');
delete(F)

%% Comparison between FHOCP and MPC trajectories
figure(10),subplot(3,1,1),plot(tsim,xMPC(4,:),'b'),grid on, hold on, title('wx [rad/s]')
figure(10),subplot(3,1,2),plot(tsim,xMPC(5,:),'b'),grid on, hold on, title('wy [rad/s]')
figure(10),subplot(3,1,3),plot(tsim,xMPC(6,:),'b'),grid on, hold on, title('wz [rad/s]')

figure(10),subplot(3,1,1),plot(tsim,xsim(4,:),'k'),grid on, hold on, title('wx [rad/s]')
figure(10),subplot(3,1,2),plot(tsim,xsim(5,:),'k'),grid on, hold on, title('wy [rad/s]')
figure(10),subplot(3,1,3),plot(tsim,xsim(6,:),'k'),grid on, hold on, title('wz [rad/s]')

figure(10),subplot(3,1,1),plot(tsim,x_ref(4,1:end-N_MPC),'r--'),grid on, hold on, title('wx [rad/s]')
figure(10),subplot(3,1,2),plot(tsim,x_ref(5,1:end-N_MPC),'r--'),grid on, hold on, title('wy [rad/s]')
figure(10),subplot(3,1,3),plot(tsim,x_ref(6,1:end-N_MPC),'r--'),grid on, hold on, title('wz [rad/s]')

legend('MPC','FHOCP','REF')

figure(11),subplot(4,1,1),plot(tsim,xMPC(1,:),'b'),grid on, hold on, title('q1')
figure(11),subplot(4,1,2),plot(tsim,xMPC(2,:),'b'),grid on, hold on, title('q2')
figure(11),subplot(4,1,3),plot(tsim,xMPC(3,:),'b'),grid on, hold on, title('q3')
figure(11),subplot(4,1,4),plot(tsim,q4_mpc(1,:),'b'),grid on, hold on, title('q4')

figure(11),subplot(4,1,1),plot(tsim,xsim(1,:),'k'),grid on, hold on, title('q1')
figure(11),subplot(4,1,2),plot(tsim,xsim(2,:),'k'),grid on, hold on, title('q2')
figure(11),subplot(4,1,3),plot(tsim,xsim(3,:),'k'),grid on, hold on, title('q3')
figure(11),subplot(4,1,4),plot(tsim,q4_sim(1,:),'k'),grid on, hold on, title('q4')

figure(11),subplot(4,1,1),plot(tsim,x_ref(1,1:end-N_MPC),'r--'),grid on, hold on, title('q1')
figure(11),subplot(4,1,2),plot(tsim,x_ref(2,1:end-N_MPC),'r--'),grid on, hold on, title('q2')
figure(11),subplot(4,1,3),plot(tsim,x_ref(3,1:end-N_MPC),'r--'),grid on, hold on, title('q3')
figure(11),subplot(4,1,4),plot(tsim,x_ref_c(4,1:end-N_MPC),'r--'),grid on, hold on, title('q4')

legend('MPC','FHOCP','REF')

% Plot control effort
figure(666),subplot(3,1,1),plot(tsim,uMPC(:,1)),grid on, hold on, title('Input acceleration [m/s^2] with MPC vs FHOCP')
figure(666),subplot(3,1,2),plot(tsim,uMPC(:,2)),grid on, hold on
figure(666),subplot(3,1,3),plot(tsim,uMPC(:,3)),grid on, hold on

legend('MPC','FHOCP')