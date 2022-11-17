function [xdot] = adcmodel(t,x,u,d,param)

Ic  = param.chaser_inertia ;
It  = param.target_inertia ;
a   = param.a ;
mu  = param.mu ;

q     =  x(1:4)   ;
wc    =  x(5:7)   ;

alfa  =  x(8:10)  ;
wt    =  x(11:13) ;

% Normalize quaternion vector
q = q/norm(q) ;

% Build skewsymmetric matric of wc
OM = [0, wc(3), -wc(2), wc(1);
      -wc(3), 0, wc(1), wc(2);
      wc(2), -wc(1), 0, wc(3);
      -wc(1), -wc(2), -wc(3), 0] ;

% Integrate Euler equations for the target in LVLH frame
itx = It(1, 1) ; ity = It(2, 2) ; itz = It(3, 3) ;

N = sqrt(mu/a^3) ;

ddalfax = - N / itx * (itz - itx - ity) * wt(2) - N^2 / itx * (itz - ity) * alfa(1) ;
ddalfay = - N / ity * (itx + ity - itz) * wt(1) - N^2 / ity * (itz - itx) * alfa(2) ;
ddalfaz = 0 ;

ddalfa = [ddalfax; ddalfay; ddalfaz] ;

% Integrate Euler equations for the chaser in body frame
q_dot = 0.5 * OM * q ;

icx = Ic(1, 1) ; icy = Ic(2, 2) ; icz = Ic(3, 3) ;

dom_x = ((icy - icz) / icx) * wc(2) * wc(3) + u(1) / icx + d(1) ;
dom_y = ((icz - icx) / icy) * wc(1) * wc(3) + u(2) / icy + d(2) ;
dom_z = ((icx - icy) / icz) * wc(2) * wc(1) + u(3) / icz + d(3) ;

dom = [dom_x; dom_y; dom_z] ;

xdot = [q_dot; dom; wt; ddalfa]  ;

end

