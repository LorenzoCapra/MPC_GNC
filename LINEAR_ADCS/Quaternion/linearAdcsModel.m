function xn = linearAdcsModel(x,u,d,param)

Ts = param.Ts;
Jinv = param.Jinv;

q = x(1:4)/norm(x(1:4));
w = x(5:7);

wdot = Jinv*u + eye(3)*d;

qdot = zeros(4,1);
qdot(1:3,1) = -0.5 * cross(w, q(1:3)) + 0.5 * q(4) * w;
qdot(4,1) = -0.5 * w' * q(1:3);

xn(5:7) = w + wdot*Ts;
xn(1:4) = q + qdot*Ts;
