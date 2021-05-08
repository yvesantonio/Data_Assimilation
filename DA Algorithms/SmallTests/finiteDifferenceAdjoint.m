function [ J_M ] = finiteDifferenceAdjoint(t, x, M)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

n = length(x);
I = eye(n);

J_M = zeros(n,n);
delta = 1e-1;

for i=1:n
    x1 = x + I(:,i)*delta;
    
    fx = M(t, x);
    fx1 = M(t, x1);
    
    diff = fx1 - fx;
    J_M(:, i) = diff/delta;    
end
% J_M = J_M';
end

