clear;
close all;

cd ../dev

% parameters
sigma   = 10;
rho     = 28;
beta    = 8/3;

% % % equation
% % f = @(t,x) [	-sigma*(x(1) - x(2));
% %                 rho*x(1) - x(2) - x(1)*x(3);
% %                 x(1)*x(2) - beta*x(3);  ];
%%% DISCRETIZATION
%   x_dot = f(x)
%   x(i+1) - x(i) = T*f(x(i))
%   x(i+1) = x(i) + T*f(x(i))
%%%
%%% ERROR, doesn't work.
%%%

%%% DISCRETIZATION as in Lawless exaple. Runge Kutta 2nd order.
M   = @(t, x, T) [  x(1) + T/2*sigma*(2*(x(2) - x(1)) + T*(rho*x(1) - x(2) - x(1)*x(3)) - T*sigma*(x(2) - x(1)));
                    x(2) + T/2*(rho*x(1) - x(2) - x(1)*x(3) + rho*(x(1) + T*sigma*(x(2) - x(1))) - x(2) - T*(rho*x(1) - x(2) - x(1)*x(3)) - (x(1) + T*sigma*(x(2) - x(1)))*(x(3) + T*(x(1)*x(2) - beta*x(3))));
                    x(3) + T/2*(x(1)*x(2) - beta*x(3) + (x(1) + T*sigma*(x(2) - x(1)))*(x(2) + T*(rho*x(1) - x(2) - x(1)*x(3))) - beta*(x(3) + T*(x(1)*x(2) - beta*x(3))));
                 ];
             
% % % % % % % % % % % % % % % % % % J_M is difficult, I generated it with Maple.
J_M = @(t, x, T) [	1 + T^2*sigma^2/2 + sigma/2*T*(-2 + T*(rho - x(3)))                                                                                             -(sigma*T/2*(T*sigma + T - 2))                                                                          -T^2/2*sigma*x(1);
                    -1/2*(T*(sigma*(beta*x(3) - 2*x(1)*x(2) + x(2)^2)*T^2 + ((rho - x(3))*sigma - (beta + 1)*x(3) + 2*x(1)*x(2) + rho)*T - 2*rho + 2*x(3)))         1 + sigma/2*(beta*x(3) + x(1)^2 - 2*x(1)*x(2))*T^3 + ((rho - x(3))*sigma - x(1)^2 + 1)*T^2/2 - T        -1/2*(T*(beta*sigma*(x(1) - x(2))*T^2 + ((-beta - sigma - 1) *x(1) + x(2)*sigma)*T + 2*x(1)));
                    -T*(sigma*((-rho/2 + x(3)/2 - 1/2)*x(2) + x(1)*(rho - x(3)))*T^2 + (x(2)*sigma/2 + x(2) + (x(3) - rho)*x(1))*T - x(2))                          1/2*((sigma*((rho - x(3) + 1)*x(1) - 2*x(2))*T^2 + ((-x(1) + 2*x(2))*sigma - 2*x(1))*T + 2*x(1))*T)     1 + (x(1)*sigma*(x(1) - x(2))*T^3)/2 + (-x(1)^2 + beta)*T^2/2 - T*beta; ];

H   = @(t, x, r) [  x(1) + random('normal', 0, r);
                    x(2) + random('normal', 0, r);
                    x(3) + random('normal', 0, r); ];
            
J_H = @(t, x) eye(3);

%%% RUN the model, to get the true result (given the discretization).
%   Use the T = 0.05;

% real
T = 0.05;
FinalTime = 5;
K = FinalTime/T +1;
x = ones(3, K);

x_0b = [ 1.2; 1.2; 1.2; ];
x_b = x_0b.*eye(3)*ones(3, K);

for i = 2:K
    x(:,i) = M(i,x(:,i-1), T);
    x_b(:,i) =	M(i,x_b(:,i-1), T);
end

%%%%% TESTS: 
% Set the options for the optimization algorythms
max_it  = 10;
tol = 1e-6;

%%% Test 1: assimilation every step.
FinalAssimilationTime = 2;
K_ass = FinalAssimilationTime/T;
y_1 = zeros(3, K_ass);

r   = 1;
% P_0b = zeros(3,3); % The book says to neglect the first term of J
P_0b = 4*(x_0b-x(:,1)).^2.*eye(3); % But the software asks for it. Invented formula.
R    = r*eye(3);    % Assume this small covariance matrix.

for i = 1:K_ass
    y_1(:,i) = H(i, x(:,i), r);
end

% % % % options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'SpecifyObjectiveGradient', true, 'MaxIterations', max_it, 'OptimalityTolerance', tol);
% % % % [ x_0a_0, Jx_0a_0] = da_var_4DVarStrongSimple(x_0b, y_1, @(t,x) M(t,x,T), @(t,x) J_M(t,x,T), @(t,x) H(t, x, 0), J_H, P_0b, R, @(f, x0) fminunc(f, x0, options));

[ x_0a_1, Jx_0a_1] = da_var_4DVarStrongSimple(x_0b, y_1, @(t,x) M(t,x,T), @(t,x) J_M(t,x,T), @(t,x) H(t, x, 0), J_H, P_0b, R, @(f,x_0) FletcherReeves(f, x_0, tol, max_it));
% Observations in every time instant.
[ x_0a_2, Jx_0a_2] = da_var_4DVarStrong(x_0b, y_1, @(t,x) M(t,x,T), @(t,x) J_M(t,x,T), @(t,x) H(t, x, 0), J_H, P_0b, R, @(f,x_0) FletcherReeves(f, x_0, tol, max_it)); %, 1:K_ass);

% x_a_1, x_a_1_bis MUST give the same result.
x_a_1 = ones(3, K);
x_a_1(:,1) = x_0a_1;

x_a_2 = ones(3, K);
x_a_2(:,1) = x_0a_2;

for i = 2:K
    x_a_1(:,i) = M(i,x_a_1(:,i-1), T);
    x_a_2(:,i) = M(i,x_a_2(:,i-1), T);
end

t = 0:K-1;
ty = 0:K_ass-1;

figure
subplot(3,1,1);
plot(t, x(1,:), 'r', ty, y_1(1,:), '*g', t, x_b(1,:), 'k.', t, x_a_1(1,:), 'b', t, x_a_2(1,:), 'b--');
subplot(3,1,2);
plot(t, x(2,:), 'r', ty, y_1(2,:), '*g', t, x_b(2,:), 'k.', t, x_a_1(2,:), 'b', t, x_a_2(2,:), 'b--');
subplot(3,1,3);
plot(t, x(3,:), 'r', ty, y_1(3,:), '*g', t, x_b(3,:), 'k.', t, x_a_1(3,:), 'b', t, x_a_2(3,:), 'b--');
legend('real', 'measurements', 'background', '4DVar', '4DVar')


%%% Test 2: assimilation every 2 (or 4  etc.) time steps, like the book suggests.
T = 0.05;
FinalAssimilationTime = 2;
K_ass = FinalAssimilationTime/T;

const = 2;
ty  = 1:const:K_ass;
y_2 = zeros(3,length(ty));

for i = ty    
    y_2(:,(i+const-1)/const) = y_1(:,i);
end

[ x_0a_3, Jx_0a_3] = da_var_4DVarStrong(x_0b, y_2, @(t,x) M(t,x,T), @(t,x) J_M(t,x,T), @(t,x) H(t, x, 0), J_H, P_0b, R, @(f,x_0) FletcherReeves(f, x_0, tol, max_it), ty);
x_a_3 = ones(3, K);
x_a_3(:,1) = x_0a_3;
for i = 2:K
    x_a_3(:,i) = M(i,x_a_3(:,i-1), T);
end
ty = ty-1;
figure
subplot(3,1,1);
plot(t, x(1,:), 'r', ty, y_2(1,:), '*g', t, x_b(1,:), 'k.', t, x_a_3(1,:), 'b');
subplot(3,1,2);
plot(t, x(2,:), 'r', ty, y_2(2,:), '*g', t, x_b(2,:), 'k.', t, x_a_3(2,:), 'b');
subplot(3,1,3);
plot(t, x(3,:), 'r', ty, y_2(3,:), '*g', t, x_b(3,:), 'k.', t, x_a_3(3,:), 'b');
cd ../test