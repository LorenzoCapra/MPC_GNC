function [xdot] = trackmodel(t,x,u,d,param)

n = param.n ;

% State space matrices
A = [0, 0, 0, 1, 0, 0;
     0, 0, 0, 0, 1, 0;
     0, 0, 0, 0, 0, 1;
     3*n^2, 0, 0, 0, 2*n, 0;
     0, 0, 0, -2*n, 0, 0;
     0, 0, -n^2, 0, 0, 0] ;
 
Bu = [0, 0, 0; 
      0, 0, 0;
      0, 0, 0;
      1, 0, 0;
      0, 1, 0;
      0, 0, 1] ;

Bd  = [0, 0, 0; 
       0, 0, 0;
       0, 0, 0;
       1, 0, 0;
       0, 1, 0;
       0, 0, 1] ;

% C = [1, 1, 1, 0, 0, 0] ;

% D = 0 ;

% d = [randn(1); randn(1); randn(1)] * 1e-4;

xdot = A*x+Bu*u+Bd*d;