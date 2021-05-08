function [ x0_ensemble ] = ensembleInit(x0b, P0b, m)
% ensembleInit generates an ensemble of points starting from the current 
% expectation  and the related uncertainty
%   x0b nx1 the current expectation
%   P0b nxn the covariance of the current expectation
%   m the dimension of the wanted ensemble

% INITIALIZATION REFERRED TO:
% Stochastic Methods for Sequential Data Assimilation in Strongly Nonlinear Systems
% Pham, 2000

n = length(x0b);

if(m<n)
    x0_ensemble = mvnrnd(x0b, P0b, m)';
    return
end

L = chol(P0b);
r1 = rank(L'); % r1 is called r' in the article, while r+1 is what I call m.
% Remembering r - r2 (>)= r1
% m - 1 - r2 = r1
% m - 1 - r1 = r2
r2 = m - r1 - 1;

C = [random('normal', 0, 1, r2,r2); zeros(m-r2, r2)];
% C = [eye(r2); zeros(m-r2, r2)];
Gamma = [ones(m, 1) C];

% r = m - 1 = r1 + r2;
z = 1; % Pick +1 or -1 with probability 1/2...
Omega = z;
% k from 2 to r - r2
% r - r2 = r1
for k = 2:r1
    z = random('normal', 0, 1, k, 1);
    % z must lie on the unitary sphere.
    z = z(:)/norm(z);
    s = sign(z(k));
    z1 = [z(1:end-1); z(end) + s];
    
    % %     H = null(z1.');

    H = eye(k) - 1/(1+z(end)*s)*(z1*z1');
    H = H(:,1:end-1);        
    
    Omega = [H*Omega z];
end

% % % Omega
% % % Omega'*Omega

H = eye(m);
for i = 0:r2
    c = Gamma(:,1);
    Gamma = Gamma(:,2:end);
    z = c/norm(c);
    s = sign(z(end));
    z1 = [z(1:end-1); z(end) + s];
    
% %     Hi = null(z1.');
    
    Hi = eye(m-i) - 1/(1+z(end)*s)*(z1*z1');
    Hi = Hi(:,1:end-1);
    
    Gamma = Hi'*Gamma;
    H = H*Hi;
end

y = sqrt(m)*H*Omega;
x0_ensemble = x0b + L'*y';
end