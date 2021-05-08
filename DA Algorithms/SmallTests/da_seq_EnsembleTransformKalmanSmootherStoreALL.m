function [ x_a_mean, x_a, x_a_meanL ] = da_seq_EnsembleTransformKalmanSmootherStoreALL(x_en, y, M, H, R, L, t_obs, varargin)
% The most simple Ensemble Kalman Filter (deterministic, also known as
% Unscent Kalman Filter).
%   x_en n*l  matrix ensemble time 0
%   y   measurements time 1 to K
%   M   f(t, x) model
%   H   f(t, x) observer
%   Q n*n covariance matrix (noise on the model)
%   R m*m covariance matrix (noise on the observer)
%   L 
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
infl = @notInflate;
if(nargin > 7)
        if(strcmp(varargin{1}, 'Inflation'))
            lambda = varargin{2};
            if(lambda > 1)
                infl = @inflate;
            end
        end
end

x_f = zeros(n,l);
x_a = zeros(n,l,K);
x_a_mean = zeros(n,K);
y_pred  = zeros(m,l);

x_a_meanL = cell(1,L);
for i=1:L
    x_a_meanL{i} = zeros(n,K);
end

x_a(:,:,1) = x_en;
x_a_mean(:,1) = mean(x_en,2);

R12 = R^-1/2;
U = eye(l);
I = eye(l);

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
        for h=1:l            
            y_pred(:,h) = H(i, x_f(:,h));
        end        
        x_f_mean = mean(x_f, 2);
        y_pred_mean = mean(y_pred, 2);
        
        X_anomaly = sq1*(x_f - x_f_mean);
        Y_anomaly = sq1*(y_pred - y_pred_mean);
        inno = y(:,j) - y_pred_mean;
        S = R12*Y_anomaly;
        T1 = I + S'*S;
        T = pinv(T1);
        
        w_a = T1\(S'*(R12*inno));
        
        x_a(:,:,i) = x_f_mean + X_anomaly*(w_a + sq*(T^(1/2)*U));
        x_a_mean(:,i) = mean(x_a(:,:,i),2);        
        x_a(:,:,i) = infl(x_a(:,:, i), x_a_mean(:,i), lambda);
        
        for ll = 1:L
            if(i-ll > 0)
                X_a_previous_anomaly = 1/sq*(x_a(:,:,i-ll) - x_a_mean(:,i-ll));
                x_a(:,:,i-ll) = x_a_mean(:,i-ll) + X_a_previous_anomaly*(w_a + sq*(T^(1/2)*U));
                x_a_mean(:,i-ll) = mean(x_a(:,:,i-ll),2);
                x_a_meanL{ll}(:,i-ll) = x_a_mean(:,i-ll);                
            end
        end
        
        j = j + 1;
    else
        x_a(:,:, i) = x_f;
        x_a_mean(:,i) = mean(x_a(:,:,i),2);
    end
end
    
%    Store lasts
for h=1:L
    x_a_meanL{h}(:, K) = x_a_mean(:, K);
end

for ll = 1:L-1
    for h=ll+1:L
        x_a_meanL{h}(:, K-ll) = x_a_meanL{ll}(:, K-ll);
    end
end

end

function [ x_a_infl ] = inflate(x_a, x_a_mean, lambda)
    x_a_infl = x_a_mean + lambda * (x_a - x_a_mean);
end

function [ x_a_infl ] = notInflate(x_a, x_a_mean, lambda)
    x_a_infl = x_a;    
end