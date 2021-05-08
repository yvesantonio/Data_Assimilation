clear;
close all;

M   = [ 1       -0.4;
        0.1     0.3; ];
Model   = @(i,x) M*x;

H   = [ 1	1];
Observer   = @(i,x) H*x;

x_b = [ 12, 175 ]';
x_0 = [ 10, 142 ]';

% % % % %%% This is a test of soundness. With this (P_0b, R) and very small tol combo we
% % % % %%% MUST get x_0_real, whateber x_b.
% % % % %%% WARNING: set R different from 0 when calling the algorythm.
P_0b = zeros(2,2);
R = 0;

P_0b   = [  5	0;
            0   50; ];
R   = 11;

x(:,1)      =   x_0;
for(i=2:3)
    x(:,i) = Model(i, x(:,i-1));
end

y(1)   = Observer(1,x(:,1)) + random('normal', 0, R);
y(2)   = Observer(2,x(:,2)) + random('normal', 0, R);
y(3)   = Observer(3,x(:,3)) + random('normal', 0, R);

% % % % % % % % % % % % % % R   = 1;
[ x_a, J_xa ] = da_var_4DVarStrongSimple(x_b, y, Model, @(k,xk) M, Observer, @(k,xk) H, P_0b, R, 1e-6, 100);

[ x_a_2, J_xa_2 ] = da_var_4DVarStrongDiscrete(x_b, y, Model, @(k,xk) M, Observer, @(k,xk) H, P_0b, R, 1e-6, 100);

[ x_a_3, J_xa_3 ] = da_var_4DVarStrongDiscrete(x_b, y(1:2:3), Model, @(k,xk) M, Observer, @(k,xk) H, P_0b, R, 1e-6, 100, [1,3]);

x_1(:,1) = x_a;
for i=2:3
    x_1(:,i) = Model(i, x_1(:,i-1));
end
figure
subplot(2,1,1);
plot(1:3, x(1,:), 1:3, x_1(1,:))
subplot(2,1,2);
plot(1:3, x(2,:), 1:3, x_1(2,:))


x_3(:,1) = x_a_3;
for i=2:3
    x_3(:,i) = Model(i, x_3(:,i-1));
end
figure
subplot(2,1,1);
plot(1:3, x(1,:), 1:3, x_3(1,:))
subplot(2,1,2);
plot(1:3, x(2,:), 1:3, x_3(2,:))