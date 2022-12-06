function xdot = AdcsLinearModel(t,x0,u,d,param)

Jinv = param.Jinv;
ws   = x0(5:7);

A   = zeros(3,3) ;
 
B   = Jinv;

Bd  = eye(3);

q = x0(1:4)/norm(x0(1:4));
% Build skewsymmetric matric of ws
OM = [0, ws(3), -ws(2), ws(1);
      -ws(3), 0, ws(1), ws(2);
       ws(2), -ws(1), 0, ws(3);
      -ws(1), -ws(2), -ws(3), 0] ;

% Integrate Euler equations for the chaser in body frame
q_dot = 0.5 * OM * q ;

wdot = A*ws + B*u + Bd*d;

xdot = [q_dot; wdot];