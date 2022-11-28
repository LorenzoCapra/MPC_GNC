function [J,xpred,ypred]=fulltrack_cost_fun(x0,upred,d,N,Q,R,x_ref,ind_sim,Ts,param)

upred       =   [upred;upred(end-1,:)];
xpred       =   zeros(size(x0,1),N+1);
ypred       =   zeros(size(x0,1),N+1);

xpred(:,1)  =   x0;
ypred(:,1)  =   x0;

slack       =   upred(end,1);

% State and output predictions
J           =   (ypred(:,1)-x_ref(:,ind_sim))'*Q*(ypred(:,1)-x_ref(:,ind_sim));

for ind_pred = 2:N
    u                   =   upred(ind_pred-1,:);
    [~, xdot]           =   ode45(@(t,x)fulltrackmodel(t,x,u,d,param), [0 Ts], xpred(:, ind_pred-1)) ;
    xpred(:,ind_pred)   =   xdot(end,:)' ;
    ypred(:,ind_pred)   =   xpred(:,ind_pred);
    J                   =   J + (ypred(:,ind_pred)-x_ref(:,ind_sim+ind_pred-1))'*Q*(ypred(:,ind_pred)-x_ref(:,ind_sim+ind_pred-1))+...
                            upred(ind_pred-1,:)*R*upred(ind_pred-1,:)';
end

J=J+1e6*slack;