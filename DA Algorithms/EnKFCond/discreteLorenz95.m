function [ x_next ] = discreteLorenz95(t, x, F, T)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

T2 = T/2;

k1 = T*Lorenz95(t, x ,F);
k2 = T*Lorenz95(t+T2, x + k1/2, F);
k3 = T*Lorenz95(t+T2, x + k2/2, F);
k4 = T*Lorenz95(t+T, x + k3, F);

x_next = x + 1/6 * (k1 + 2*k2 + 2*k3 + k4);
end

