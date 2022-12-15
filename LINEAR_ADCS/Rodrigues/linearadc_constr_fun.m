function [Constr,Constr_eq,xpred]=linearadc_constr_fun(x0,upred,d,N,umax,param)

wmax = param.wmax;

Constr_eq                                   =   [];
Constr                                      =   zeros(2*size(upred,1)+2*N+1,6);
Constr(1:size(upred,1),1:3)                 =   upred-umax;
Constr(size(upred,1)+1:2*size(upred,1),1:3) =   -upred-umax;

upred       =   [upred;upred(end-1,:)];
slack       =   upred(end,1);
xpred       =   zeros(size(x0,1),size(upred,1));
xpred(:,1)  =   x0;

Constr(1:size(upred,1),4:6)                 =   (xpred(1:3,:)-wmax)' ;
Constr(size(upred,1)+1:2*size(upred,1),4:6) =   (-xpred(1:3,:)-wmax)';

% State and output predictions

for ind_pred = 2:N+1
    u                   =   upred(ind_pred-1,:)';
    xn                  =   linearAdcsModel(xpred(:,ind_pred-1),u,d,param);
    xpred(:,ind_pred)   =   xn ;
end

Constr(end,1)   =   -slack;