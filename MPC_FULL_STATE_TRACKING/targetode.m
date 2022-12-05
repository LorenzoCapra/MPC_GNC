function [xdot] = targetode(t,x,param)

It  = param.target_inertia ;

qt      =   x(1:4);
wt      =   x(5:7);

%% Rotation
% Normalize quaternion vector
qt = qt/norm(qt) ;

% Build skewsymmetric matric of wc
OM = [0, wt(3), -wt(2), wt(1);
      -wt(3), 0, wt(1), wt(2);
      wt(2), -wt(1), 0, wt(3);
      -wt(1), -wt(2), -wt(3), 0] ;

% Integrate Euler equations for the chaser in body frame
qt_dot = 0.5 * OM * qt ;

itx = It(1, 1) ; ity = It(2, 2) ; itz = It(3, 3) ;

dom_x = ((ity - itz) / itx) * wt(2) * wt(3);
dom_y = ((itz - itx) / ity) * wt(1) * wt(3);
dom_z = ((itx - ity) / itz) * wt(2) * wt(1);

dom = [dom_x; dom_y; dom_z] ;

%% Return state derivative
xdot = [qt_dot; dom];
