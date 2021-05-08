% Lorenz Equations Script
clear; close all;

% parameters
sigma   = 10;
rho     = 28;
beta    = 8/3;

% equation
x_dot = @(t,x)  [   -sigma*(x(1) - x(2));
                    rho*x(1) - x(2) - x(1)*x(3);
                    x(1)*x(2) - beta*x(3);  ];

% solution
T   = 100;
x0  = [1; 1; 1;];

options = odeset('RelTol', 1e-8, 'MaxStep', 0.05);
[t, x]  = ode23(x_dot, [0,T], x0, options);

figure
plot3(x(:,1), x(:,2), x(:,3));
xlabel('x');
ylabel('y');
zlabel('z');

figure
plot(t, x(:,1));
figure
plot(t, x(:,2));
figure
plot(t, x(:,3));

derivatives = zeros(3,length(t));
for i=1:length(t)
    derivatives(:,i) = x_dot(0,x(i,:));
end
figure
plot(t, derivatives (1,:));
hold on
plot(t, derivatives (2,:));
hold on
plot(t, derivatives (3,:));