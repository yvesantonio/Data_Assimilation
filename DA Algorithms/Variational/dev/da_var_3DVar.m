function [ x_a, J_xa] = da_var_3DVar(x_b, y, B, R, H, f_min)
% 3DVar gets x_a where x_a is optimal according to the 3D-Var cost function:
%   J   =   1/2*(x - x_b)'* B^-1 *(x - x_b) + 1/2*(y - H*x)'*R^-1*(y - H*x);
%   DJ  =   B^-1*(x - x_b) + H'*R^-1*(y - H*x)
%   READ da_var4DStrong documentation

B1  =   pinv(B);
R1  =   pinv(R);

[x_a, J_xa] = f_min(@(x) func(x, x_b, y, B1, R1, H), x_b)
end

function [J, DJ] = func(x, x_b, y, B1, R1, H)
J   =   1/2*(x - x_b)'* B1 *(x - x_b) + 1/2*(y - H*x)'*R1*(y - H*x);
DJ  =   B1*(x - x_b) - H'*R1*(y - H*x);
end