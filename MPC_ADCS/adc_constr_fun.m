function [Constr,Constr_eq,xpred]=adc_constr_fun(x0,upred,d,N,umax,Ts,param)

wmax = param.wmax;

Constr_eq                                   =   [];
Constr                                      =   zeros(2*size(upred,1)+2*N+1,6);
Constr(1:size(upred,1),1:3)                 =   upred-umax;
Constr(size(upred,1)+1:2*size(upred,1),1:3) =   -upred-umax;

upred       =   [upred;upred(end-1,:)];
slack       =   upred(end,1);
xpred       =   zeros(size(x0,1),size(upred,1));
xpred(:,1)  =   x0;

%xpred_mpc   =   zeros(size(x0,1),20+1);
Constr(1:size(upred,1),4:6)                 =   (xpred(5:7,:)-wmax)' ;
Constr(size(upred,1)+1:2*size(upred,1),4:6) =   (-xpred(5:7,:)-wmax)';

% State and output predictions

for ind_pred = 2:N+1
    u                   =   upred(ind_pred-1,:);
    [~, xdot]           =   ode45(@(t,x)adcmodel(t,x,u,d,param), [0 Ts], xpred(:, ind_pred-1)) ;
    xpred(:,ind_pred)   =   xdot(end,:)' ;
    %xpred(:,ind_pred)  =   Ad*xpred(:,ind_pred-1)+Bd*u';
end

Constr(end,1)   =   -slack;
