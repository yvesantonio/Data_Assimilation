function [ x_0a, Jx_0a] = da_var_4DVarStrongSimple(x_0b, y, M, J_M, H, J_H, P_0b, R, f_min)
% 4DVar gets x_a where x_a is optimal according to the strong constrained 4D-Var cost function:
%   J   =   1/2*(x - x_b)'* P^-1 *(x - x_b) + 1/2*sum((y(k) - H(k,x(k)))'*R^-1*(y(k) - H(k,x(k))));
%   x_b (nx1)   background of x_0,
%   y   (mxK)   measurements vector
%   M(k,xk)     data generation model:  x(k+1)  = M(k, xk);
%   J_M(k,xk)	Jacbian of M; TODO: make this optional.
%   H(k,xk)     data observation model: y(k)    = H(k, xk);
%   J_H(k,xk)	Jacbian of H; TODO: make this optional.
%
%   f_min       Optimization function: [x_opt, f_opt] = @(J, x0) f_min(J,
%   x0)
%   where J is such thut [Jx, DJx] = J(x)
%   Tolerances, max iterations and other options have to be set outside this function.
%
%   Please, set 'SpecifyObjectiveGradient' to true when using Matlab
%   Optimization toolbex (such as fminunc).

K = length(y(1,:));
n = length(x_0b);
m = length(y(:,1));
if (length(size(R)) < 3)
    s	=   size(R);
    R1	=  zeros(s(1), s(2), K);
    for i = 1:K
        R1(:,:,i) = pinv(R);
    end
else
    R1 = R;
    for i=1:K
        R1(:,:,i) = pinv(R(:,:,i));
    end
end

P1  =   pinv(P_0b);
[x_0a, Jx_0a, counter, error, xks, fks] = f_min(@(x_0) func(x_0, x_0b, n, K, m, y, M, J_M, H, J_H, P1, R1), x_0b);
end

function [J, DJ] = func(x_0, x_0b, n, K, m, y, M, J_M, H, J_H, P1, R1)
t = 1:K;
% Compute x with forward model

x = zeros(n, K);
x(:,1) = x_0;
for i = 2:K
    x(:,i) = M(i, x(:,i-1));
end

% Compute cost function J
to_sum  = zeros(1,K);
d       = zeros(m,K);
for i = t
    d(:,i) = y(:,i) - H(i, x(:,i));
    to_sum(i) = d(:,i)'*R1(:,:,i)*d(:,i);
end
J = 1/2*(x_0 - x_0b)'* P1 *(x_0 - x_0b) + 1/2*sum(to_sum);

% Compute final condition pK for adjoint model
p       =   zeros(n,K);
p(:,K)  =   -J_H(K, x(K))'*R1(:,:,K)*d(:,K);

% Compute p with adjoint model
for i = K-1:-1:2
    p(:,i) = -J_H(i, x(:,i))'*R1(:,:,i)*d(:,i) + J_M(i, x(:,i-1))'*p(:,i+1);
end
DJ	=	P1*(x(:,1)-x_0b) + J_M(1, x(:,1))'*p(:,2) -J_H(1, x(:,1))'*R1(:,:,1)*d(:,1);
end