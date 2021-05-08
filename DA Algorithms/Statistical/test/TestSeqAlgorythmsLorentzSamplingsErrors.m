clear;
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
FinalAssimilationTime = 20;
K = FinalAssimilationTime/T +1;
x = ones(3, K);

%%%x_0b = H(i, x(:,1), 4);
x_0b = [1.2; 1.2; 1.2;];
x_b = x_0b.*eye(3)*ones(3, K);

for i = 2:K
    x(:,i) = M(i,x(:,i-1), T);
    x_b(:,i) =	M(i,x_b(:,i-1), T);
end

y_sets = 1001;
y = zeros(3, K-1, y_sets);

r   = 4;
P_0b = r*eye(3);
R    = r*eye(3);    % Assume this covariance matrix.
Q    = zeros(3,3);
parfor i = 2:K
    for j=1:y_sets
        y(:, i-1, j) = H(i, x(:,i), r);
    end
end

max_samples = 10;
c = 1:max_samples;

x_0b_en = ensembleInit(x_0b, P_0b, 50);

ARMSE_EKF = zeros(y_sets, max_samples);
ARMSE_EnKF = zeros(y_sets, max_samples);
ARMSE_EnKFInfl = zeros(y_sets, max_samples);
ARMSE_ETKF = zeros(y_sets, max_samples);

RMSE_EKF = zeros(y_sets, max_samples);
RMSE_EnKF = zeros(y_sets, max_samples);
RMSE_EnKFInfl = zeros(y_sets, max_samples);
RMSE_ETKF = zeros(y_sets, max_samples);

for i=c
    h = 1:i:K;
    x_cur = x(:,1:h(end));
    parfor j=1:y_sets
        y_cur = y(:,h(1:end-1),j);
                
        [ x_a_EKF, P_a_EKF] = da_seq_ExtendedKalmanFilter(x_0b, y_cur, @(t,x) M(t,x,T), @(t,x) J_M(t,x,T), @(t,x) H(t, x, 0), J_H, P_0b, Q, R, h);
        [ ARMSE_EKF(j, i), RMSE] = averageRootMeanSquareError(x_a_EKF, x_cur);
        RMSE_EKF(j, i) = RMSE(end);
        
        [ x_a_EnKF, x_a_EnKF_en] = da_seq_StochasticEnsembleKalmanFilter(x_0b_en, y_cur, @(t,x) M(t,x,T), @(t,x) H(t, x, 0), R, h);
        [ ARMSE_EnKF(j, i), RMSE] = averageRootMeanSquareError(x_a_EnKF, x_cur);
        RMSE_EnKF(j,i) = RMSE(end);
        
        [ x_a_EnKFInfl, x_a_EnKFInfl_en] = da_seq_StochasticEnsembleKalmanFilter(x_0b_en, y_cur, @(t,x) M(t,x,T), @(t,x) H(t, x, 0), R, h, 'Inflation', 1.1);
        [ ARMSE_EnKFInfl(j, i), RMSE] = averageRootMeanSquareError(x_a_EnKFInfl, x_cur);
        RMSE_EnKFInfl(j,i) = RMSE(end);

        [ x_a_ETKF, x_a_ETKF_en] = da_seq_EnsembleTransformKalmanFilter(x_0b_en, y_cur, @(t,x) M(t,x,T), @(t,x) H(t, x, 0), R, h);
        [ ARMSE_ETKF(j, i), RMSE] = averageRootMeanSquareError(x_a_ETKF, x_cur);
        RMSE_ETKF(j, i) = RMSE(end);
    end
    i
    save samplingErrors
end

ARMSE_EKF = (ARMSE_EKF);
ARMSE_EnKF = mean(ARMSE_EnKF);
ARMSE_EnKFInfl = mean(ARMSE_EnKFInfl);
ARMSE_ETKF = mean(ARMSE_ETKF);

RMSE_EKF = mean(RMSE_EKF);
RMSE_EnKF = mean(RMSE_EnKF);
RMSE_EnKFInfl = mean(RMSE_EnKFInfl);
RMSE_ETKF = mean(RMSE_ETKF);

figure
plot(c, ARMSE_EKF(c), 'g*', c, ARMSE_EnKFInfl(c), 'b*-', c, ARMSE_ETKF(c), 'c*-');
legend('EKF ARMSE', 'EnKF ARMSE Infl', 'ETKF ARMSE');

figure
plot(c, ARMSE_EnKF(c), 'b*--' , c, ARMSE_EnKFInfl(c), 'b*-');
legend('EnKF ARMSE No Inflation', 'EnKF ARMSE Optimal Inflation');

figure
plot(c, RMSE_EKF, 'g*', c, RMSE_EnKF(c), 'b*--', c, RMSE_EnKFInfl(c), 'b*-', c, RMSE_ETKF(c), 'c*-');
legend('EKF last RMSE', 'EnKF ARMSE No Infl', 'EnKF  last RMSE Optimal Infl', 'ETKF last RMSE');

cd ../test
save samplingErrors
