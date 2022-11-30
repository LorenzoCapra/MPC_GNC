function [xdot] = targetode(t,x,param)

n   = param.n ;
It  = param.target_inertia ;

alfa    =   x(1:3);
wt      =   x(4:6);

%% Rotation
% Integrate Euler equations for the target in LVLH frame
itx = It(1, 1) ; ity = It(2, 2) ; itz = It(3, 3) ;

ddalfax = - n / itx * (itz - itx - ity) * wt(2) - n^2 / itx * (itz - ity) * alfa(1) ;
ddalfay = - n / ity * (itx + ity - itz) * wt(1) - n^2 / ity * (itz - itx) * alfa(2) ;
ddalfaz = 0 ;

ddalfa = [ddalfax; ddalfay; ddalfaz] ;

%% Return state derivative
xdot = [wt;ddalfa];
