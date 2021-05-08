function [ x_a_mean, x_a ] = da_seq_bundleIterativeEnsembleKalmanSmoother(x_en, y, M, H, R, options)
% The Iterative Ensemble Kalman Smoother (Filter as a particular case).
%   x_en n*l  matrix ensemble time 0
%   y   measurements time 1 to K
%   M   f(t, x) model
%   H   f(t, x) observer
%   R m*m covariance matrix (noise on the observer)
%   options: generate with the relative function

n = length(x_en(:,1));
m = length(y(:,1));
l = length(x_en(1,:));
K = length(y(1,:)) + 1;

if(nargin > 5)
    [ L, S, lambda, epsilon, fMin, beta, selector, loc ] = readOptions(options);
else
    [ L, S, lambda, epsilon, fMin, beta, selector, loc ] = readOptions(bIEnKSOptions());
end

if(lambda <= 1)
    infl = @notInflate;
else
    infl = @inflate;
end

x_f = zeros(n,l,K);
x_a = zeros(n,l,K);
x_a_mean = zeros(n,K);

x_a(:,:,1) = x_en;
x_a_mean(:,1) = mean(x_en,2);

U = eye(l);

sq = sqrt(l-1);
epsilon = epsilon*sq;
sq1 = 1/sq;

y = [y, zeros(m, L)];

if(loc == 0)
    analysis = @globalAnalysis;
    R12 = R^-(1/2);
else
    analysis = @localAnalysis;
    sel = cell(1, n);
	R12 = cell(1, n);
    for h=1:n
        s = selector(h);
        sel{h} = s;        
        R12{h} = (s*R*s')^-(1/2);
    end
    selector = sel;
end

i_start = 2+L;
if(beta == 0)
    objective = @objectiveSDA;
    
    % Preprocessing: NOT TO WASTE ANY MEASUREMENT.
    % SAME WINDOW SIZE, all measurement are part of J(w)
    % False for to make the code foldable.
    for i = 2
        % Step 1: prediction
        % Project state ahead
        for h=1:l
            x_f(:,h,i) = M(i, x_a(:,h,i-1));
        end
        % Step 2: measurements update
        x_f_mean = mean(x_f(:,:,i), 2);
        X_anomaly = sq1*(x_f(:,:,i) - x_f_mean);        
        
        x_a(:,:, i) = analysis(objective, fMin, x_f_mean, X_anomaly, L, L+1, M, H, R12, y, epsilon, n, m, l, L+2, sq, U, beta, selector);
        x_a_mean(:,i) = mean(x_a(:,:,i),2);
        
        % % %
        % % %     MULTIPLICATIVE INFLATION HERE
        % % %
        x_a(:,:, i) = infl(x_a(:,:, i), x_a_mean(:,i), lambda);
        
        for s=i+1:i+S-1
            for h=1:l
                x_a(:,h, s) = M(i, x_a(:,h, s-1));
            end
            x_a_mean(:,s) = mean(x_a(:,:,s),2);
        end
    end
    i_start = 2+L+S;
else
    objective = @objectiveMDA;
end

for i=i_start:S:K+L
    % Step 1: prediction
    %i
    % Project state ahead    
    for h=1:l
        x_f(:,h,i-L) = M(i, x_a(:,h,i-L-1));        
    end
    
    % Step 2: measurements update if available
    x_f_mean = mean(x_f(:,:,i-L), 2);
    X_anomaly = sq1*(x_f(:,:,i-L) - x_f_mean);

    x_a(:,:, i-L) = analysis(objective, fMin, x_f_mean, X_anomaly, L, S, M, H, R12, y, epsilon, n, m, l, i, sq, U, beta, selector);
    x_a_mean(:,i-L) = mean(x_a(:,:,i-L),2);
    
    % % %
    % % %     MULTIPLICATIVE INFLATION HERE
    % % %
    x_a(:,:, i-L) = infl(x_a(:,:, i-L), x_a_mean(:,i-L), lambda);
    
    for s=i-L+1:min(i-L+S-1, K)
        for h=1:l
            x_a(:,h, s) = M(i, x_a(:,h, s-1));
        end
        x_a_mean(:,s) = mean(x_a(:,:,s),2);
    end
end

end

function [ L, S, lambda, epsilon, fMin, beta, selector, loc ] = readOptions(options)
    L = options.L;
    S = options.S;
    
    lambda = options.Inflation;
    epsilon = options.epsilon;
    fMin = options.fMin;
    beta = options.beta;
    
    loc = options.loc;
    selector = options.selector;
end

function [ x_a_infl ] = inflate(x_a, x_a_mean, lambda)
    x_a_infl = x_a_mean + lambda * (x_a - x_a_mean);
end

function [ x_a_infl ] = notInflate(x_a, x_a_mean, lambda)
    x_a_infl = x_a;
end

function [ f, g, B ] = objectiveSDA(w, x0_mean, X0_anomaly, L, S, M, H, R12, y, epsilon, n, m, l, i, beta)
    base = L-S+2; % Matlab 1 = 0 on the book.
    k = base:L+1; % 1 2 3 ... L+1 = L on the book
    
    x_mean = zeros(n, L+1);
    yH_mean = zeros(m, L+1);
    
    x0 = x0_mean + X0_anomaly*w;    
    
    x_mean(:, 1) = x0;
    yH_mean(:, 1) = H(i, x_mean(:,1));
    
    for kk = 2:L+1
        x_mean(:, kk) = M(i+kk, x_mean(:, kk-1));
        yH_mean(:, kk) = H(i+kk, x_mean(:, kk));
    end
    
    f = 0.5*(w'*w);
    for j=k
        tmp = (y(:, j) - yH_mean(:, j));
        SMatrix = R12*tmp;
        f = f + 0.5*(SMatrix'*SMatrix);
    end
    
    if(nargout > 1)
        x = zeros(n, l, L);
        x(:,:, 1) = x0 + epsilon*X0_anomaly;
        
        yH = zeros(m, l, L);
        yH_mean = zeros(m, L);
        Y_anomaly = zeros(m, l, L);
        
        for h=1:l
            yH(:,h, 1) = H(i+kk, x(:,h, 1));
        end
        yH_mean(:,1) = mean(yH(:,:, 1), 2);
        
        for kk = 2:L+1
            for h=1:l
                x(:, h, kk) = M(i+kk, x(:, h, kk-1));
                yH(:, h, kk) = H(i+kk, x(:, h, kk));
            end
            yH_mean(:, kk) = mean(yH(:,:,kk), 2);
        end
        
        g = w;
        for j=k
            Y_anomaly(:,:, j) = 1/epsilon*(yH(:,:, j) - yH_mean(:, j));
            SMatrix = R12*Y_anomaly(:,:, j);
            g = g - (SMatrix'*(R12*(y(:, j) - yH_mean(:, j))));
        end
        
        if(nargout > 2)
            B = eye(l);
            for j=k
                SMatrix = R12*Y_anomaly(:,:, j);
                B = B + SMatrix'*SMatrix;
            end
        end
    end
end

function [ f, g, B ] = objectiveMDA(w, x0_mean, X0_anomaly, L, S, M, H, R12, y, epsilon, n, m, l, i, beta)
    k = 1:L+1; % 1 2 3 ... L+1 = L on the book
    
    x_mean = zeros(n, L+1);
    yH_mean = zeros(m, L+1);
    
    x0 = x0_mean + X0_anomaly*w;    
    
    x_mean(:, 1) = x0;
    yH_mean(:, 1) = H(i, x_mean(:,1));
    
    for kk = 2:L+1
        x_mean(:, kk) = M(i+kk, x_mean(:, kk-1));
        yH_mean(:, kk) = H(i+kk, x_mean(:, kk));
    end
    
    f = 0.5*(w'*w);
    for j=k
        tmp = (y(:, j) - yH_mean(:, j));
	SMatrix = R12*tmp;
        f = f + beta(j)*0.5*(SMatrix'*SMatrix);
    end
    
    if(nargout > 1)
        x = zeros(n, l, L);
        x(:,:, 1) = x0 + epsilon*X0_anomaly;
        
        yH = zeros(m, l, L);
        yH_mean = zeros(m, L);
        Y_anomaly = zeros(m, l, L);
        
        for h=1:l
            yH(:,h, 1) = H(i+kk, x(:,h, 1));
        end
        yH_mean(:,1) = mean(yH(:,:, 1), 2);
        
        for kk = 2:L+1
            for h=1:l
                x(:, h, kk) = M(i+kk, x(:, h, kk-1));
                yH(:, h, kk) = H(i+kk, x(:, h, kk));
            end
            yH_mean(:, kk) = mean(yH(:,:,kk), 2);
        end
        
        g = w;
        for j=k
            Y_anomaly(:,:, j) = 1/epsilon*(yH(:,:, j) - yH_mean(:, j));
            SMatrix = R12*Y_anomaly(:,:, j);
            g = g - beta(j)*(SMatrix'*(R12*(y(:, j) - yH_mean(:, j))));
        end
        
        if(nargout > 2)
            B = eye(l);
            for j=k
                SMatrix = R12*Y_anomaly(:,:, j);
                B = B + beta(j)*(SMatrix'*SMatrix);
            end
        end
    end
end

function [ x_a ] = globalAnalysis(objective, fMin, x_f_mean, X_anomaly, L, S, M, H, R12, y, epsilon, n, m, l, i, sq, U, beta, selector)
w = zeros(l,1);
w_a = fMin(@(w) objective(w, x_f_mean, X_anomaly, L, S, M, H, R12, y(:, i-L-1:i-1), epsilon, n, m, l, i, beta), w);
[tmp, tmp1, B] = objective(w_a, x_f_mean, X_anomaly, L, S, M, H, R12, y(:, i-L-1:i-1), epsilon, n, m, l, i, beta);

x_a =  x_f_mean + X_anomaly*w_a + sq*(X_anomaly*(pinv(B)^(1/2)*U));
end

function [ x_a ] = localAnalysis(objective, fMin, x_f_mean, X_anomaly, L, S, M, H, R12, y, epsilon, n, m, l, i, sq, U, beta, selector)
x_a = zeros(n, l);
parfor h=1:n
    w = zeros(l,1);
    sel = selector{h};
    y_cur = sel*y(:, i-L-1:i-1);
    H_cur = @(t,x) sel*H(t,x);
    m_cur = length(y_cur(:,1));
    w_a = fMin(@(w) objective(w, x_f_mean, X_anomaly, L, S, M, H_cur, R12{h}, y_cur, epsilon, n, m_cur, l, i, beta), w);
    [tmp, tmp1, B] = objective(w_a, x_f_mean, X_anomaly, L, S, M, H_cur, R12{h}, y_cur, epsilon, n, m_cur, l, i, beta);
    
    x_a(h,:) =  x_f_mean(h) + X_anomaly(h,:)*w_a + sq*(X_anomaly(h,:)*(pinv(B)^(1/2)*U));
end
end
