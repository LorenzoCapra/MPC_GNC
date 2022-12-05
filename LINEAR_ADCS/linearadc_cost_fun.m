function [J,xpred,ypred]=linearadc_cost_fun(x0,upred,d,N,Q,R,x_ref,ind_sim,A,B)

upred       =   [upred;upred(end-1,:)];
xpred       =   zeros(size(x0,1),N+1);
ypred       =   zeros(size(x0,1),N+1);

xpred(:,1)  =   x0;
ypred(:,1)  =   x0;

slack       =   upred(end,1);

% State and output predictions
if ind_sim == 1
    ind_sim = ind_sim + 1;
end
J           =   (ypred(:,1)-x_ref(:,ind_sim-1))'*Q*(ypred(:,1)-x_ref(:,ind_sim-1));

for ind_pred = 2:N
    u                   =   upred(ind_pred-1,:);
    xpred(:,ind_pred)   =   A*xpred(:,ind_pred-1)+B*u';
    ypred(:,ind_pred)   =   xpred(:,ind_pred);
    J                   =   J + (ypred(:,ind_pred)-x_ref(:,ind_sim+ind_pred-1))'*Q*(ypred(:,ind_pred)-x_ref(:,ind_sim+ind_pred-1))+...
                            upred(ind_pred-1,:)*R*upred(ind_pred-1,:)';
end

J=J+1e6*slack;