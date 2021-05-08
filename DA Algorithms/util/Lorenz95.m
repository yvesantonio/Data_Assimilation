function [ x_dot ] = Lorenz95(t, x, F)
% The Lorenz-95 model
% t hasn't got a role in he computations, it is a legacy choise.
% x nx1 current set of parameters
% F energy injection

x_dot = zeros(length(x),1);
x_dot(1) = (x(2) - x(end-1))*x(end) - x(1) + F;
x_dot(2) = (x(3) - x(end))*x(1) - x(2) + F;
x_dot(3:end-1) = (x(4:end) - x(1:end-3)).*x(2:end-2) - x(3:end-1) + F;
x_dot(end) = (x(1) - x(end-2))*x(end-1) - x(end) + F;

end

