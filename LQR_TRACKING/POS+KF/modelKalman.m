function xdot = modelKalman(t,x,u,y,sys)

Ak       = sys.A;
Bk       = sys.B;

xdot = Ak*x + Bk*[u;y];