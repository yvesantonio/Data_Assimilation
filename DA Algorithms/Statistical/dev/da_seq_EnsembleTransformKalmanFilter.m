function [ x_a_mean, x_a ] = da_seq_EnsembleTransformKalmanFilter(x_en, y, M, H, R, t_obs, varargin)
% The most simple Ensemble Kalman Filter (deterministic, also known as
% Unscent Kalman Filter).
%   x_en n*l  matrix ensemble time 0
%   y   measurements time 1 to K
%   M   f(t, x) model
%   H   f(t, x) observer
%   Q n*n covariance matrix (noise on the model)
%   R m*m covariance matrix (noise on the observer)
%   Options: {'Inflation', lambda}, {'Localization', selector(i)}

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
R1 = R^-1;
analysis = @globalAnalysis;
selector = 0;
if(nargin > 6)
    for i=1:2:length(varargin)
        if(strcmp(varargin{i}, 'Inflation'))
            lambda = varargin{i+1};
            if(lambda > 1)
                infl = @inflate;                
            end
        elseif(strcmp(varargin{i}, 'Localization'))
            analysis = @localAnalysis;
            sel = varargin{i+1};
            selector = cell(1, n);
            R1 = cell(1, n);
            for h=1:n
                s = sel(h);
                selector{h} = s;
                R1{h} = (s*R*s')^-1;
            end
        end
    end
end

x_f = zeros(n,l);
x_a = zeros(n,l,K);
x_a_mean = zeros(n,K);
y_pred  = zeros(m,l);

x_a(:,:,1) = x_en;
x_a_mean(:,1) = mean(x_en,2);

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
        X_anomaly = sq1*(x_f - x_f_mean);
        
        x_a(:,:,i) = analysis(x_f_mean, X_anomaly, R1, y(:, j), y_pred, n, sq, sq1, I, U, selector);
        
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

function [ x_a ] = globalAnalysis(x_f_mean, X_anomaly, R1, y, y_pred, n, sq, sq1, I, U, selector)
y_pred_mean = mean(y_pred, 2);

Y_anomaly = sq1*(y_pred - y_pred_mean);
inno = y - y_pred_mean;
            
T = (I + Y_anomaly'*R1*Y_anomaly)^-1;

w_a = T*Y_anomaly'*R1*inno;

x_a = x_f_mean + X_anomaly*(w_a + sq*T^(1/2)*U);
end

function [ x_a ] = localAnalysis(x_f_mean, X_anomaly, R1, y, y_pred, n, sq, sq1, I, U, selector)
parfor h=1:n
    sel = selector{h};
    y_local = sel*y_pred;
    
    y_local_mean = mean(y_local, 2);
    
    Y_anomaly = sq1*(y_local - y_local_mean);
    inno = sel*y - y_local_mean;
    T = (I +  Y_anomaly'*R1{h}*Y_anomaly)^-1;
    w_a = T*Y_anomaly'*R1{h}*inno;
    
    x_a(h,:) = x_f_mean(h) + X_anomaly(h,:)*(w_a + sq*T^(1/2)*U);
end
end