function [J,xpred,ypred]=track_cost_fun(x0,upred,N,Q,R,Ad,Bd,Cd,Dd,x_ref,ind_sim)

upred       =   [upred;upred(end-1,:)];
xpred       =   zeros(6,N+1);
ypred       =   zeros(6,N+1);
xpred(:,1)  =   x0;
ypred(:,1)  =   Cd*x0+Dd*upred(1,:)';
slack       =   upred(end,1);

% State and output predictions
J           =   (ypred(:,1)-x_ref(:,ind_sim))'*Q*(ypred(:,1)-x_ref(:,ind_sim));

for ind_pred = 2:N%+1
    u                   =   upred(ind_pred-1,:);
    xpred(:,ind_pred)   =   Ad*xpred(:,ind_pred-1)+Bd*u';
    ypred(:,ind_pred)   =   Cd*xpred(:,ind_pred)+Dd*upred(ind_pred,:)';
    J                   =   J + (ypred(:,ind_pred)-x_ref(:,ind_sim+ind_pred-1))'*Q*(ypred(:,ind_pred)-x_ref(:,ind_sim+ind_pred-1))+...
                            upred(ind_pred-1,:)*R*upred(ind_pred-1,:)';
end

J=J+1e6*slack;