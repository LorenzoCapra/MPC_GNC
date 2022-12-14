function xdot = model(t,x,u,d,param)

n = param.n ;

% State space matrices for realtive motion
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

xdot = A*x+Bu*u+Bd*d;