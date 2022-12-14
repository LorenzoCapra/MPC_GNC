function [xF,yF,P] = KalmanFilter(sys,Y,u,W,V,P0,x0)

%%%%% NOT WORKING PROPERLY %%%%%

% INPUTS:
% - sys : linear state space system
% - Y   : measurements/output signal
% - u   : input signal
% - W   : variance of the process noise w
% - V   : variance of measurements noise v
% - P0  : initial condition variance
% - x0  : initial state vector
%
% OUTPUT:
% - xF      : prediction of the state
% - yF      : predicted output
% - P       : variance of the state error
%
% Author: Lorenzo Capra - Politecnico di Milano %

A = sys.A;
B = sys.B;
C = sys.C;

% Prediction step
xp  = A*x0 + B*u;
Pp  = A*P0*A' + W;

KF  = Pp*C' / (C*Pp*C' + V);

% Compute error
yp  = C*xp ;
e   = Y - yp;

% Correction step
x   = xp + KF*e ;
P   = Pp - KF*C*Pp;

% Return estimated state and output
xF  = x;
yF  = C*x;
    
