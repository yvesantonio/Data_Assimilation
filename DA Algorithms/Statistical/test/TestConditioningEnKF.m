hclear;
close all;

cd ../dev

% parameters
sigma	= 10;
rho     = 28;
beta    = 8/3;

M   = @(t, x, T) [  x(1) + T/2*sigma*(2*(x(2) - x(1)) + T*(rho*x(1) - x(2) - x(1)*x(3)) - T*sigma*(x(2) - x(1)));
                    x(2) + T/2*(rho*x(1) - x(2) - x(1)*x(3) + rho*(x(1) + T*sigma*(x(2) - x(1))) - x(2) - T*(rho*x(1) - x(2) - x(1)*x(3)) - (x(1) + T*sigma*(x(2) - x(1)))*(x(3) + T*(x(1)*x(2) - beta*x(3))));
                    x(3) + T/2*(x(1)*x(2) - beta*x(3) + (x(1) + T*sigma*(x(2) - x(1)))*(x(2) + T*(rho*x(1) - x(2) - x(1)*x(3))) - beta*(x(3) + T*(x(1)*x(2) - beta*x(3))));
                 ];

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
FinalAssimilationTime = 2;
K = FinalAssimilationTime/T +1;
x = ones(3, K);

%%%x_0b = H(i, x(:,1), 4);
x_0b = [1.2; 1.2; 1.2;];
x_b = x_0b.*eye(3)*ones(3, K);

for i = 2:K
    x(:,i) = M(i,x(:,i-1), T);
    x_b(:,i) =	M(i,x_b(:,i-1), T);
end

%%% Test 1: assimilation every step.
y_sets = 101;
y = zeros(3, K-1, y_sets);

r   = 4;
P_0b = r*eye(3);
R    = r*eye(3);    % Assume this small covariance matrix.
Q    = zeros(3,3);
parfor i = 2:K
    for j=1:y_sets
        y(:, i-1, j) = H(i, x(:,i), r);        
    end
end

min = 4;
max = 50;

RConditioningYY = zeros(y_sets, max, K);
ConditioningYY = zeros(y_sets, max, K);

% % % % % % % plot(1:K, RMSE_b, 'k*-', 1:K, RMSE_EKF, 'g*-')
% 3 states, start with an ensemble of 4 points and go up to 100
ens = min:max;
for i=ens
    x_0b_en = ensembleInit(x_0b, P_0b, i);
    parfor j=1:y_sets        
        [ x_a_mean, x_a_en, ConditioningYY(j, i, :), RConditioningYY(j, i, :) ] = da_seq_StochasticEnKFCOND(x_0b_en, y(:,:,j), @(t,x) M(t, x, T) , @(t,x) H(t, x, 0), R);
    end
end

% Plot condition number ongoing (pick a random measurements set), few
% elements & lot of elements
figure
plot(1:K, squeeze(RConditioningYY(55, 4, :)), 1:K, squeeze(RConditioningYY(55, 50, :)));
legend('rcond, 4 elements ensemble', 'rcond, 50 elements ensemble');

% Plot condition number ongoing (pick a random measurements set), few
% elements & lot of elements
figure
plot(1:K, squeeze(RConditioningYY(5, 4, :)), 1:K, squeeze(RConditioningYY(5, 50, :)));
legend('rcond, 4 elements ensemble', 'rcond, 50 elements ensemble');

% Plot condition number ongoing (pick a random measurements set), few
% elemens & lot of elements
figure
plot(1:K, squeeze(RConditioningYY(100, 4, :)), 1:K, squeeze(RConditioningYY(100, 50, :)));
legend('rcond, 4 elements ensemble', 'rcond, 50 elements ensemble');

RCYY_average = mean(RConditioningYY, 3);

figure
plot(ens, RCYY_average(55, ens), ens, RCYY_average(5, ens), ens, RCYY_average(100, ens));

RCYY_ongoing = mean(RCYY_average);

figure
plot(ens, RCYY_ongoing(ens));


CYY_max = zeros(y_sets, max);
clear max
parfor i=ens
	l = squeeze(ConditioningYY(:,i,:));
    CYY_max(:, i) = max(l');
end

figure
plot(ens, CYY_max(55, ens), ens, CYY_max(5, ens), ens, CYY_max(100, ens));

CYY_average = mean(CYY_max);

figure
plot(ens, CYY_average(ens));

cd ../test