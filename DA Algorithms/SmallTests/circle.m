function [ x ] = circle( i, x, tetha_dot, T)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
A = [   0 tetha_dot;
        -tetha_dot 0 ];

k1 = A*x;
k2 = A*(x + 1/2*k1*T);
k3 = A*(x + 1/2*k2*T);
k4 = A*(x + k3*T);

x = x + T/6*(k1 + 2*k2 + 2*k3 + k4);
end

