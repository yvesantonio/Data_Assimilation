clear;
close all;

cd ../dev

%%%%
%%% TRY TO REPLICATE THE SAME EXACT ALGORYTHM FOR THE STATISTICAL DA
%%%%

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
             
% % % % % % % % % % % % % % % % % % NOTE: J_M is difficult and long to compute
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

%%% Test 1: assimilation every step.
FinalAssimilationTime = 2;
K_ass = FinalAssimilationTime/T;
y_1 = zeros(3, K_ass-1);

r   = 2;
P_0b = 4*(x_0b-x(:,1)).^2.*eye(3);
R    = r*eye(3);    % Assume this small covariance matrix.
Q    = zeros(3,3);
for i = 2:K_ass
    y_1(:,i-1) = H(i, x(:,i), r);    
end

x_0b_en = ensembleInit(x_0b, P_0b, 50);

tic
[ x_a_MLETKF10, x_a_MLETKF10_en] = da_hyb_MaximumLikelihoodEnsembleTransformKalmanFilter(x_0b_en, y_1, @(t,x) M(t,x,T), @(t,x) H(t, x, 0), R, 1:K_ass, 'MaxIter', 5);
toc

% USELESS: perfectly identical.
% % % % tic
% % % % [ x_a_MLETKF, x_a_MLETKF_en] = da_hyb_MaximumLikelihoodEnsembleTransformKalmanFilter(x_0b_en, y_1, @(t,x) M(t,x,T), @(t,x) H(t, x, 0), R);
% % % % toc

x_a_MLETKF10 = [x_a_MLETKF10 zeros(3,K-K_ass)];
% % % % x_a_MLETKF = [x_a_MLETKF zeros(3,K-K_ass)];

for i = K_ass+1:K
    x_a_MLETKF10(:,i) = M(i,x_a_MLETKF10(:,i-1), T);
    %%%%x_a_MLETKF(:,i) = M(i,x_a_MLETKF(:,i-1), T);
end

t = 0:K-1;
ty = 1:K_ass-1;
figure
subplot(3,1,1);
plot(t, x(1,:), 'r', ty, y_1(1,:), '*g', t, x_b(1,:), 'k.', t, x_a_MLETKF10(1,:), 'g'); %, t, x_a_MLETKF(1,:), 'g--');
subplot(3,1,2);
plot(t, x(2,:), 'r', ty, y_1(2,:), '*g', t, x_b(2,:), 'k.', t, x_a_MLETKF10(2,:), 'g'); %, t, x_a_MLETKF(2,:), 'g--');
subplot(3,1,3);
plot(t, x(3,:), 'r', ty, y_1(3,:), '*g', t, x_b(3,:), 'k.', t, x_a_MLETKF10(3,:), 'g'); %, t, x_a_MLETKF(3,:), 'g--');
legend('real', 'measurements', 'background', 'MLETKF 5 iterations/cycle', 'MLETKF');

figure
tt = 0:K_ass-1;
subplot(3,1,1);
plot(tt, x(1,tt+1), 'r', ty, y_1(1,:), '*g', tt, x_b(1,tt+1), 'k.', tt, x_a_MLETKF10(1,tt+1), 'g'); %, tt, x_a_MLETKF(1,tt+1), 'g--');
subplot(3,1,2);
plot(tt, x(2,tt+1), 'r', ty, y_1(2,:), '*g', tt, x_b(2,tt+1), 'k.', tt, x_a_MLETKF10(2,tt+1), 'g'); %, tt, x_a_MLETKF(2,tt+1), 'g--');
subplot(3,1,3);
plot(tt, x(3,tt+1), 'r', ty, y_1(3,:), '*g', tt, x_b(3,tt+1), 'k.', tt, x_a_MLETKF10(3,tt+1), 'g'); %, tt, x_a_MLETKF(3,tt+1), 'g--');
legend('real', 'measurements', 'background', 'MLETKF 5 iterations/cycle', 'MLETKF');

%%% Test 2: assimilation every 2 (or 4  etc.) time steps, like the book suggests.
T = 0.05;
FinalAssimilationTime = 2;
K_ass = FinalAssimilationTime/T;

const = 4;
ty  = 1:const:K_ass;
y_2 = zeros(3,length(ty));

for i = ty
    y_2(:,(i+const-1)/const) = y_1(:,i);
end

tic
[ x_a_MLETKF10_2, x_a_MLETKF10_2_en] = da_hyb_MaximumLikelihoodEnsembleTransformKalmanFilter(x_0b_en, y_2, @(t,x) M(t,x,T), @(t,x) H(t, x, 0), R, ty, 'MaxIter', 5);
toc

K_ass = length(x_a_MLETKF10_2(1,:));
x_a_MLETKF10_2 = [x_a_MLETKF10_2 zeros(3,K-K_ass)];

for i = K_ass:K
    x_a_MLETKF10_2(:,i) = M(i,x_a_MLETKF10_2(:,i-1), T);
end

figure
subplot(3,1,1);
plot(t, x(1,:), 'r', ty, y_2(1,:), '*g', t, x_b(1,:), 'k.', t, x_a_MLETKF10_2(1,:), 'g');
subplot(3,1,2);
plot(t, x(2,:), 'r', ty, y_2(2,:), '*g', t, x_b(2,:), 'k.', t, x_a_MLETKF10_2(2,:), 'g');
subplot(3,1,3);
plot(t, x(3,:), 'r', ty, y_2(3,:), '*g', t, x_b(3,:), 'k.', t, x_a_MLETKF10_2(3,:), 'g');
legend('real', 'measurements', 'background', 'MLETKF 5 iterations/cycle');
cd ../test