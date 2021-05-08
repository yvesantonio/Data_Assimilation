function [ x_a_mean, x_a ] = da_seq_StochasticEnsembleKalmanFilter(x_en, y, M, H, R, t_obs, varargin)
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

infl = @notInflate;
lambda = -1;
localization = 0;
selector = 0;
if(nargin > 6)
    for i=1:2:length(varargin)
        if(strcmp(varargin{i}, 'Inflation'))
            lambda = varargin{i+1};
            if(lambda > 1)
                infl = @inflate;                
            end
        elseif(strcmp(varargin{i}, 'Localization'))
            localization = 1;
            sel = varargin{i+1};
            selector = cell(1, n);
            Rl = cell(1, n);
            for h=1:n
                s = sel(h);
                selector{h} = s;
                Rl{h} = (s*R*s');
            end
        end
    end
end

x_f = zeros(n,l);
x_a = zeros(n,l,K);
x_a_mean = zeros(n,K);
u = zeros(m,l);

x_a(:,:,1) = x_en;
x_a_mean(:,1) = mean(x_en,2);

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
        x_f_mean = mean(x_f, 2);
        X_anomaly = sq1*(x_f - x_f_mean);
        if(~localization)
            y_set   = zeros(m,l);
            y_pred  = zeros(m,l);
            for h=1:l
                u(:,h) = mvnrnd(zeros(m,1), R);
                y_set(:,h) = y(:,j) + u(:,h);
                y_pred(:,h) = H(i, x_f(:,h));
            end
            u_mean = mean(u, 2);
            y_pred_mean = mean(y_pred,2);
            
            Y_anomaly = sq1*(y_pred - u - y_pred_mean + u_mean);
            % Compute the Kalman gain
            K_g = X_anomaly*Y_anomaly'*pinv(Y_anomaly*Y_anomaly');
            
            % Update estimate with measurements
            x_a(:,:,i) = x_f + K_g*(y_set - y_pred);
        else
            y_cur = y(:,j);  
            yHcur = zeros(m,l);
            parfor hl=1:l
                yHcur(:,hl) = H(i, x_f(:,hl));
            end
            parfor hh=1:n
                sel = selector{hh};
                y_set   = sel*zeros(m,l);
                y_pred  = sel*zeros(m,l);
                u = sel*zeros(m,l);
                z = sel*zeros(m,1);
                for h=1:l
                    u(:,h) = mvnrnd(z, Rl{hh});
                    y_set(:,h) = sel*y_cur + u(:,h);
                    y_pred(:,h) = sel*yHcur(:,h);
                end
                u_mean = mean(u, 2);
                y_pred_mean = mean(y_pred,2);
                
                Y_anomaly = sq1*(y_pred - u - y_pred_mean + u_mean);
                
                K_g = X_anomaly*Y_anomaly'*pinv(Y_anomaly*Y_anomaly');
                K_g = K_g(hh,:);
                % Update estimate with measurements
                x_a(hh,:,i) = x_f(hh,:) + K_g*(y_set - y_pred);
            end
        end
        j = j + 1;
    else
        x_a(:,:, i) = x_f;
    end
    x_a_mean(:,i) = mean(x_a(:,:,i),2);
    
    % % %
    % % %     MULTIPLICATIVE INFLATION HERE
    % % %
    x_a(:,:, i) = infl(x_a(:,:, i), x_a_mean(:,i), lambda);    
end
end

function [ x_a_infl ] = inflate(x_a, x_a_mean, lambda)
    x_a_infl = x_a_mean + lambda * (x_a - x_a_mean);
end

function [ x_a_infl ] = notInflate(x_a, x_a_mean, lambda)
    x_a_infl = x_a;    
end