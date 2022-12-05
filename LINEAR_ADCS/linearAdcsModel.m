function [xdot] = linearAdcsModel(t,x,u,d,param)

Jinv = param.Jinv;

A   = [0.5*eye(3,3) ,  zeros(3,3);
       zeros(3,3)   ,  zeros(3,3)] ;
 
B   = [zeros(3,3); Jinv];

Bd  = [zeros(3,3); eye(3)];

xdot = A*x + B*u + Bd*d;
