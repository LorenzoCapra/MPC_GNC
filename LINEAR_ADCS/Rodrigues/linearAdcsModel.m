function xn = linearAdcsModel(x,u,d,param)

Ts = param.Ts;
Jinv = param.Jinv;
J = param.chaser_inertia;

o = x(4:6);
w = x(1:3);

So = SkewSymmetric(o);
Sw = SkewSymmetric(w);

G = 0.5*(eye(3) + o*o' - So - (1+o'*o)*eye(3)/2);

odot = G*w;
wdot = Jinv*Sw*J*w + Jinv*u + Jinv*d;

xn(1:3) = w + wdot*Ts;
xn(4:6) = o + odot*Ts;
