function [ options ] = bIEnKSOptions( varargin )
% UNTITLED Summary of this function goes here
%   Detailed explanation goes here

options.L = 5;
options.S = 2;
options.epsilon = 0.1;
options.fMin = @(f,x0) QuasiNewtonBFGS(f, x0, 1e-6, 100);
options.Inflation = 1;
options.beta = 0;
options.selector = 0;
options.loc = 0;

i = 1;
while i < length(varargin)
    opt = varargin{i};
    if(strcmp(opt, 'L'))
        options.L = varargin{i+1};
        i = i+2;
    elseif(strcmp(opt, 'S'))
        options.S = varargin{i+1};
        i = i+2;
    elseif(strcmp(opt, 'epsilon'))
        options.epsilon = varargin{i+1};
        i = i+2;
	elseif(strcmp(opt, 'OptAlg'))
        options.fMin = varargin{i+1};
        i = i+2;
    elseif(strcmp(opt, 'Inflation'))
        options.Inflation = varargin{i+1};
        i = i+2;
    elseif(strcmp(opt, 'Filter'))
        options.L = 1;
        options.S = 1;
    elseif(strcmp(opt, 'MDA'))
        options.beta = varargin{i+1};
        i = i+2;
	elseif(strcmp(opt, 'Localization'))
        options.loc = 1;
        options.selector = varargin{i+1};
        i = i+2;
    end
end

if(options.beta ~= 0)
    if(length(options.beta) ~= options.L+1)
        error('beta must be of length L+1')        
    elseif(norm(sum(options.beta) - options.S) > 0.1)
        error('sum(beta) must be S')
    end
end

end

