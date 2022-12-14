function xdot = modelDist(t,x,u,sys)

A       = sys.A;
B       = sys.B;

xdot = A*x + B*u;