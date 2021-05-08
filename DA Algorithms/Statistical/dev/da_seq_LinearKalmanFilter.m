function [ x_a, P_a] = da_seq_LinearKalmanFilter(x_0b, y, M, H, P_0b, Q, R, t_obs) 
% The most simple Kalman Filter, with M the model,
%   x_0b first guess for time 0
%   y m*K measurements foe time 1 to K, or for t=t_obs
%   M n*n matrix
%   H m*n matrix
%   Q n*n covariance matrix (noise on the model)
%   R m*m covariance matrix (noise on the observer)
%   t_obs vector of available sampled time starting atleast from time 1
n = length(x_0b);
m = length(y(:,1));
if(nargin < 8)
    K = length(y(1,:)) + 1;
    t_obs = 1:K;    
else
    K = t_obs(end);    
end

P_f = zeros(n,n,K);
P_a = zeros(n,n,K);
K_g = zeros(n,m,K);
x_f = zeros(n,K);
x_a = zeros(n,K);

x_a(:,1) = x_0b;
P_a(:,:,1) = P_0b;

j = 1;
for i=2:K
    % Step 1: prediction
    
    % Project state ahead
    x_f(:,i) = M*x_a(:, i-1);
    % Project the error covariance ahead
    P_f(:,:,i) = M*P_a(:,:,i-1)*M' + Q;
    
    % Step 2: measurements update IF AVAILABLE
    if i-1==t_obs(j)
        % Compute the Kalman gain
        K_g(:,:,i) = P_f(:,:,i)*H'*(H*P_f(:,:,i)*H' + R)^-1;
        % Update estimate with measurements
        x_a(:,i) = x_f(:,i) + K_g(:,:,i)*(y(:,j) - H*x_f(:,i));
        % Update the error covariance (JUST THE SIMPLEST OF THE POSSIBLE
        % FORMULAE)
        P_a(:,:,i) = (eye(n) - K_g(:,:,i)*H)*P_f(:,:,i);
        j = j+1;
    else
        x_a(:,i) = x_f(:,i);
        P_a(:,:,i) = P_f(:,:,i);
    end
end

end