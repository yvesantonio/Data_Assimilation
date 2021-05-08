function [ ARMSE, RMSE ] = averageRelativeRootMeanSquareError( x_analysis, x_real)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    n = length(x_real(:,1));
    K = length(x_real(1,:));
    v = x_analysis - x_real;
    dot_prod = zeros(1,K);
    for i=1:K
        dot_prod(i) = v(:,i)'*v(:,i);
        dot_prod(i) = dot_prod(i)/(x_real(:,i)'*x_real(:,i));
        dot_prod = dot_prod/n;
    end
    RMSE = sqrt(dot_prod);
    ARMSE = mean(RMSE);
end

