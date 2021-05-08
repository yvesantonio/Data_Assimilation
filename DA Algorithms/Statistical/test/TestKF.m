clear;
close all;

cd ../dev;
Q = 5*eye(3);
P_0b = diag([50, 4, 6]);

% Choose a stable System, not to make things difficult at this stage
M = [	0.9     -1      0.3;
        0       -1.2    -1.1;
        0.14    0.1     -0.4; ];

H = [1 1 0;
     0 0 1];
    
% % % % % % % % a = eig(M)
% % % % % % % % t = 0:0.1:2*pi;
% % % % % % % % x = cos(t);
% % % % % % % % y = sin(t);
% % % % % % % % plot(a,'*')
% % % % % % % % hold on;
% % % % % % % % plot(x, y);

K = 200;
x_0b = [150; 10; 15];
x_t  = zeros(3,K);
x_t(:,1) = [135; 8.5; 17];
y = zeros(2,K-1);

for i=2:K
    x_t(:,i) = M*x_t(:,i-1);
    y(:,i-1) = H*x_t(:,i) + [ random('normal', 0, 5); random('normal', 0, 2)];
end
R = diag([5, 2]);

[ x_a_check, P_a_check] = da_seq_LinearKalmanFilterSimple(x_0b, y, M, H, P_0b, Q, R);
[ x_a, P_a] = da_seq_LinearKalmanFilter(x_0b, y, M, H, P_0b, Q, R);

t = 0:K-1;
figure
subplot(3,1,1);
plot(t, x_t(1,:), 'r', t, x_a(1,:), t, x_a_check(1,:), 'g')
subplot(3,1,2);
plot(t, x_t(2,:), 'r', t, x_a(2,:), t, x_a_check(2,:), 'g')
subplot(3,1,3);
plot(t, x_t(3,:), 'r', t, x_a(3,:), t, x_a_check(3,:), 'g')

P_0b
P_a_end = P_a(:,:,K)
eig(P_a_end)
P_a_end = P_a_check(:,:,K)
eig(P_a_end)

t_2 = 1:2:K;
y_2 = y(:,t_2);

[ x_a_2, P_a_2] = da_seq_LinearKalmanFilter(x_0b, y_2, M, H, P_0b, Q, R, t_2);

t = 0:K-1;
figure
subplot(3,1,1);
plot(t, x_t(1,:), 'r', t(1:K-1), x_a_2(1,:))
subplot(3,1,2);
plot(t, x_t(2,:), 'r', t(1:K-1), x_a_2(2,:))
subplot(3,1,3);
plot(t, x_t(3,:), 'r', t(1:K-1), x_a_2(3,:))
P_a_end = P_a_2(:,:,K-1)
eig(P_a_end)

t_3 = 2:5:K;
y_3 = y(:,t_3);

[ x_a_3, P_a_3] = da_seq_LinearKalmanFilter(x_0b, y_3, M, H, P_0b, Q, R, t_3);

t = 0:K-1;
figure
subplot(3,1,1);
plot(t, x_t(1,:), 'r', t(1:K-3), x_a_3(1,:))
subplot(3,1,2);
plot(t, x_t(2,:), 'r', t(1:K-3), x_a_3(2,:))
subplot(3,1,3);
plot(t, x_t(3,:), 'r', t(1:K-3), x_a_3(3,:))
P_a_end = P_a_3(:,:,K-3)
eig(P_a_end)

cd ../test