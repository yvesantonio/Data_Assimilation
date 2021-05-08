function [ x_a_mean, x_a, condYY ] = da_seq_StochasticEnKFCOND(x_en, y, M, H, R, t_obs, varargin)
% The most simple Ensemble Kalman Filter
%   x_en n*l  matrix ensemble time 0
%   y   measurements time 1 to K
%   M f(t, x) model
%   H f(t, x) observer 
%   Q n*n covariance matrix (noise on the model)
%   R m*m covariance matrix (noise on the observer)
%   Options: {'Inflation', lambda}

n = length(x_en(:,1));
m = length(y(:,1));
l = length(x_en(1,:));

if(nargin < 6)
    K = length(y(1,:)) + 1;
    t_obs = 1:K;    
else
    K = t_obs(end);
end

lambda = -1;
if(nargin > 6)
    for i=1:2:length(varargin)
        if(strcmp(varargin{i}, 'Inflation'))
            lambda = varargin{i+1};
        end
    end
end

condYY = zeros(1, K);
x_f = zeros(n,l);
x_a = zeros(n,l,K);
x_a_mean = zeros(n,K);
y_set   = zeros(m,l);
y_pred  = zeros(m,l);
u = zeros(m,l);

x_a(:,:,1) = x_en;
x_a_mean(:,1) = mean(x_en,2);

sq = sqrt(l-1);
sq1 = 1/sq;

I = eye(m);

j = 1;
for i=2:K
    % Step 1: prediction
    
    % Project state ahead
    for h=1:l
        x_f(:,h) = M(i, x_a(:,h,i-1));
    end
    
    % Step 2: measurements update if available
    if i-1==t_obs(j)        
        for h=1:l
            u(:,h) = mvnrnd(zeros(m,1), R)';
            y_set(:,h) = y(:,j) + u(:,h);
            y_pred(:,h) = H(i, x_f(:,h));
        end
        x_f_mean = mean(x_f, 2);
        u_mean = mean(u, 2);
        y_pred_mean = mean(y_pred,2);
        
        X_anomaly = sq1*(x_f - x_f_mean);
        Y_anomaly = sq1*(y_pred - y_pred_mean - u + u_mean);
        Matrix = Y_anomaly*Y_anomaly';
        try
            condYY(i) = cond(Matrix, 2);
        catch
            condYY(i) = NaN;
        end
        Matrix = pinv(Matrix);
        % Compute the Kalman gain
        K_g = (X_anomaly*Y_anomaly')*Matrix;
        
        % Update estimate with measurements
        for h=1:l            
            x_a(:,h,i) = x_f(:,h) + K_g*(y_set(:,h) - y_pred(:,h));
        end
        j = j + 1;
    else
        x_a(:,:, i) = x_f;
    end
    x_a_mean(:,i) = mean(x_a(:,:,i),2);
    
    % % %
    % % %     MULTIPLICATIVE INFLATION HERE
    % % %
    
    if(lambda > 1)
        x_a(:,:, i) = x_a_mean(:,i) + lambda * (x_a(:,:, i) - x_a_mean(:,i));
    end
end
end