function [ x_next ] = discreteLorenz95(t, x, F, T)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[t, x] = ode45(@(t, y) Lorenz95(t,y,F), [0 T], x);
x_next = x(end,:)';
end

