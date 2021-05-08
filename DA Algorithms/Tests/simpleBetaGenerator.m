function [ beta ] = simpleBetaGenerator( n )

beta = 0.5*ones(1,n);

for i=n-1:-1:2
    beta(i) = beta(i+1)/2;
end
beta(1) = beta(2);
end

