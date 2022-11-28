function [Constr,Constr_eq,xpred]=cw_constr_fun(x0,upred,N,Ad,Bd,umax)

Constr_eq                                 =   [];
Constr                                    =   zeros(2*size(upred,1)+2*N+1,3); % Input constraints;State constraints
Constr(1:size(upred,1),:)                 =   upred-umax;
Constr(size(upred,1)+1:2*size(upred,1),:) =   -upred-umax;

upred       =   [upred;upred(end-1,:)];
slack       =   upred(end,1);
xpred       =   zeros(6,N+1);
xpred(:,1)  =   x0;

% State and output predictions

for ind_pred = 2:N+1
    u                   =   upred(ind_pred-1,:);
    xpred(:,ind_pred)   =   Ad*xpred(:,ind_pred-1)+Bd*u';
end

Constr(end,1)   =   -slack;