function [xdot] = fulltrackmodel(t,x,u,d,param)

% d = [normrnd(0, 1e-5); normrnd(0, 1e-5); normrnd(0, 1e-5);
%      normrnd(0, 1e-6); normrnd(0, 1e-6); normrnd(0, 1e-6)];

n = param.n ;
Ic  = param.chaser_inertia ;

r       =   x(1:3);
v       =   x(4:6);
q       =   x(7:10);
wc      =   x(11:13);

%% Translation
xdd     =   -3*n^2 * r(1) + 2*n*v(2) + u(1) + d(1);
ydd     =   -2*n*v(1) + u(2) + d(2);
zdd     =   -n^2 * r(3) + u(3) + d(3);

rdd     =   [xdd;ydd;zdd];

%% Rotation
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

dom_x = ((icy - icz) / icx) * wc(2) * wc(3) + u(4) / icx + d(4) ;
dom_y = ((icz - icx) / icy) * wc(1) * wc(3) + u(5) / icy + d(5) ;
dom_z = ((icx - icy) / icz) * wc(2) * wc(1) + u(6) / icz + d(6) ;

dom = [dom_x; dom_y; dom_z] ;

%% Return state derivative
xdot = [v;rdd;q_dot;dom];
