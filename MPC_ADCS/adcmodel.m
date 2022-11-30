function [xdot] = adcmodel(t,x,u,d,param)

Ic  = param.chaser_inertia ;

q     =  x(1:4)   ;
wc    =  x(5:7)   ;

% Normalize quaternion vector
q = q/norm(q) ;

% Build skewsymmetric matric of wc
OM = [0, wc(3), -wc(2), wc(1);
      -wc(3), 0, wc(1), wc(2);
      wc(2), -wc(1), 0, wc(3);
      -wc(1), -wc(2), -wc(3), 0] ;

% Integrate Euler equations for the chaser in body frame
q_dot = 0.5 * OM * q ;

icx = Ic(1, 1) ; icy = Ic(2, 2) ; icz = Ic(3, 3) ;

dom_x = ((icy - icz) / icx) * wc(2) * wc(3) + u(1) / icx + d(1) ;
dom_y = ((icz - icx) / icy) * wc(1) * wc(3) + u(2) / icy + d(2) ;
dom_z = ((icx - icy) / icz) * wc(2) * wc(1) + u(3) / icz + d(3) ;

dom = [dom_x; dom_y; dom_z] ;

xdot = [q_dot; dom]  ;

end

