function [Constr,Constr_eq,xpred]=fulltrack_constr_fun(x0,upred,d,N,umax,Ts,param)

wmax = param.wmax;

Constr_eq                                   =   [];
Constr                                      =   zeros(2*size(upred,1)+2*N+1,9);
Constr(1:size(upred,1),1:6)                 =   upred-umax;
Constr(size(upred,1)+1:2*size(upred,1),1:6) =   -upred-umax;

upred       =   [upred;upred(end-1,:)];
slack       =   upred(end,1);
xpred       =   zeros(size(x0,1),size(upred,1));
xpred(:,1)  =   x0;

Constr(1:size(upred,1),7:9)                 =   (xpred(5:7,:)-wmax)' ;
Constr(size(upred,1)+1:2*size(upred,1),7:9) =   (-xpred(5:7,:)-wmax)';

% State and output predictions

for ind_pred = 2:N+1
    u                   =   upred(ind_pred-1,:);
    [~, xdot]           =   ode45(@(t,x)fulltrackmodel(t,x,u,d,param), [0 Ts], xpred(:,ind_pred-1)) ;
    xpred(:,ind_pred)   =   xdot(end,:)' ;
end

Constr(end,1)   =   -slack;