function [ J_M ] = lorentz95Adjoint(t, x, T)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

n = length(x);
I = eye(n);

J_M = zeros(n,n);

J_M(1, n-1) = -T*x(n);
J_M(1, n) = T*(x(2) - x(n-1));
J_M(1, 1) = 1 - T;
J_M(1, 2) = T*x(n);

J_M(2, n) = -T*x(1);
J_M(2, 1) = T*(x(3) - x(n));
J_M(2, 2) = 1 - T;
J_M(2, 3) = T*x(1);

for i=3:n-1    
    J_M(i, i-2) = -T*x(i-1);
    J_M(i, i-1) = T*(x(i+1) - x(i-2));
    J_M(i, i) = 1 - T;
    J_M(i, i+1) = T*x(i-1); 
end

J_M(n, n-2) = -T*x(n-1);
J_M(n, n-1) = T*(x(1) - x(n-2));
J_M(n, n) = 1 - T;
J_M(n, 1) = T*x(n-1);

end

