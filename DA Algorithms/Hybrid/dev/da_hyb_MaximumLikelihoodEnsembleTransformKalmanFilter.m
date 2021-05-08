function [ x_a_mean, x_a ] = da_hyb_MaximumLikelihoodEnsembleTransformKalmanFilter(x_en, y, M, H, R, t_obs, varargin)
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

if(nargin < 6)
    K = length(y(1,:)) + 1;
    t_obs = 1:K;    
else
    K = t_obs(end);
end

lambda = -1;
tol = 1e-8;
max_it = 100;
epsilon = 0.1;
if(nargin > 6)
    for i=1:2:length(varargin)        
        if(strcmp(varargin{i}, 'Inflation'))
            lambda = varargin{i+1};
        elseif(strcmp(varargin{i}, 'Tolerance'))
            tol = varargin{i+1};
    	elseif(strcmp(varargin{i}, 'MaxIter'))
            max_it = varargin{i+1};
        elseif(strcmp(varargin{i}, 'Epsilon'))
            epsilon = varargin{i+1};
        end
    end
end

x_f = zeros(n,l);
x_a = zeros(n,l,K);
x_a_mean = zeros(n,K);
y_pred  = zeros(m,l);

x_a(:,:,1) = x_en;
x_a_mean(:,1) = mean(x_en,2);

j = 1;

R1 = R^-1;
U = eye(l);
sq = sqrt(l-1);
for i=2:K
    % Step 1: prediction
    
    % Project state ahead
    for h=1:l
        x_f(:,h) = M(i, x_a(:,h,i-1));        
    end
    
    % Step 2: measurements update if available
    if i-1==t_obs(j)
        x_f_mean = mean(x_f, 2);
        X_anomaly = 1/sq*(x_f - x_f_mean);
        
        T = U; % WARNING: Theory: T is eye(l)
        w = zeros(l,1);
        opt_index = 0;
        
        Dw = tol*100;
        
        %%%%%%%
        %%%%%   WARNING: The optimization method is HARD CODED
        %%%%%%%
        while ( norm(Dw) >= tol || opt_index <= max_it)            
            x = x_f_mean + X_anomaly*w;
            
% % %             % Bundle
% % %             E = x + epsilon*X_anomaly;
            % Transform
            E = x + sq*X_anomaly*T;
            
            for h=1:l
                y_pred(:,h) = H(i, E(:,h));
            end
            y_pred_mean = mean(y_pred, 2);
            
% % %             % Bundle
% % %             Y_anomaly = (y_pred - y_pred_mean)/epsilon;
            % Transform
            Y_anomaly = 1/sq*(y_pred - y_pred_mean)*T^-1;
            
            inno = y(:,j) - y_pred_mean;
            grad = w - Y_anomaly'*R1*inno;
            hess = eye(l) + Y_anomaly'*R1*Y_anomaly;
            Dw = hess\grad;
            w = w - Dw;
            
            % Transform
            T = hess^(-1/2);
            opt_index = opt_index+1;
        end                     
% % %         % Bundle
% % %         T = hess^(-1/2);
        x_a(:,:, i) = x + sq*X_anomaly*T*U; 
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