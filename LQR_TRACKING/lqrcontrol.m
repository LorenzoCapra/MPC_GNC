function u = lqrcontrol(sys,x0,xref,Q,R)

% A = sys.A;
% B = sys.B;
% C = sys.C;
% D = sys.D;

% Aa = [A, zeros(6,3);
%       -C(1:3,:), zeros(3,3)];
%   
% Ba = [B; zeros(3,3)];
% 
% Ca = [C, zeros(6,3)];
% 
% sysa = ss(Aa,Ba,Ca,0);

% rank(ctrb(Aa,Ba))
% [Aa_bar, Ba_bar, Ca_bar, T, k] = ctrbf(Aa, Ba, Ca);
% n_uc = size(Aa, 1) - sum(k) % Number of uncontrollable modes is 8 - 6 = 2
% Aa_uc = Aa_bar(1:n_uc, 1:n_uc)

[K,~,~] = lqr(sys, Q, R);

u = -K*(xref-x0);