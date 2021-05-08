function [ x_a_mean, x_a ] = da_seq_LocalizedEnsembleTransformKalmanFilter(x_en, y, M, J_H, R, rho, t_obs, varargin)
% The most simple Ensemble Kalman Filter (deterministic, also known as
% Unscent Kalman Filter).
%   x_en n*l  matrix ensemble time 0
%   y   measurements time 1 to K
%   M   f(t, x) model
%   H   f(t, x) observer
%   Q n*n covariance matrix (noise on the model)
%   R m*m covariance matrix (noise on the observer)
%   Options: {'Inflation', lambda}

n = length(x_en(:,1));
m = length(y(:,1));
l = length(x_en(1,:));

if(nargin < 7)
    K = length(y(1,:)) + 1;
    t_obs = 1:K;    
else
    K = t_obs(end);
end

lambda = -1;
if(nargin > 7)
    for i=1:2:length(varargin)
        if(strcmp(varargin{i}, 'Inflation'))
            lambda = varargin{i+1};
        end
    end
end

x_f = zeros(n,l);
x_a = zeros(n,l,K);

x_a_mean = zeros(n,K);

x_a(:,:,1) = x_en;
x_a_mean(:,1) = mean(x_en,2);

R1 = R^-1;
U = eye(l);
I = eye(n);

sq = sqrt(l-1);
sq1 = 1/sq;

j = 1;
for i=2:K
    % Step 1: prediction
    
    % Project state ahead
    for h=1:l
        x_f(:,h) = M(i, x_a(:,h,i-1));        
    end
    
    % Step 2: measurements update if available
    if i-1==t_obs(j)
        y_pred = J_H*x_f; 
        
        x_f_mean = mean(x_f, 2);        
        y_pred_mean = mean(y_pred, 2);        
        
        X_anomaly = sq1*(x_f - x_f_mean);        
        
        inno = y(:,j) - y_pred_mean;
        
        B = rho.*(X_anomaly*X_anomaly');
        
        Kg = B*J_H'*pinv(J_H*B*J_H' + R);
        T = pinv(I + B*J_H'*R1*J_H);
        x_a(:,:,i) = x_f_mean + Kg*inno + sq*T^(1/2)*X_anomaly*U;
        
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