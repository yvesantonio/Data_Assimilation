function [ x_a, P_a] = da_seq_ExtendedKalmanFilterSimple(x_0b, y, M, J_M, H, J_H, P_0b, Q, R) 
% The most simple Kalman Filter, with M the model,
%   x_0b first guess time 0
%   y   measurements time 1 to K
%   M       f(t, x) model
%   J_M     f(t, x) model Jacobian
%   H m*n matrix
%   Q n*n covariance matrix (noise on the model)
%   R m*m covariance matrix (noise on the observer)

n = length(x_0b);
m = length(y(:,1));
K = length(y(1,:)) + 1;
P_f = zeros(n,n,K);
P_a = zeros(n,n,K);
K_g = zeros(n,m,K);
x_f = zeros(n,K);
x_a = zeros(n,K);

x_a(:,1) = x_0b;
P_a(:,:,1) = P_0b;

for i=2:K
    % Step 1: prediction
    
    % Project state ahead
    x_f(:,i) = M(i, x_a(:, i-1));
    % Project the error covariance ahead
    J = J_M(i, x_a(:, i-1));
    P_f(:,:,i) = J*P_a(:,:,i-1)*J' + Q;
    
    % Step 2: measurements update
    % Compute the Kalman gain
    J = J_H(i, x_a(:, i-1));
    K_g(:,:,i) = P_f(:,:,i)*J'*(J*P_f(:,:,i)*J' + R)^-1;
    % Update estimate with measurements
    x_a(:,i) = x_f(:,i) + K_g(:,:,i)*(y(:,i-1) - H(i, x_f(:,i)));
    % Update the error covariance (JUST THE SIMPLEST OF THE POSSIBLE
    % FORMULAE)
    P_a(:,:,i) = (eye(n) - K_g(:,:,i)*J)*P_f(:,:,i);
end

end