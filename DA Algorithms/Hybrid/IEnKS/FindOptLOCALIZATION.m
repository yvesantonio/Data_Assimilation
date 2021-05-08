%%%
%%% Localization: look for the optimum localization radius
%%%

%%%
%%% WARNING: Make sure to change the names of the variables and save the
%%% results of the script.
%%%

clear;
close all;

load BookTestData


%%%%%%%%%%%%%%
%%%%%%%%%%%%%% WARNING: K is changed in order to make the computation
%%%%%%%%%%%%%% faster.
%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Run 2:
LWindow = [1 5 10 15];
L_max = 10;
min = 7;
max = 15;
loc = 3:7;
%%%%%%%%%%%%%% FILE GENERATED: ETKFOptLocTest1
%%%%%%%%%%%%%%

ens = min:max;
r = loc(end);

ARMSE = zeros(r, L_max, max);
RMSE = cell(r, L_max, max);
ARelRMSE = zeros(r, L_max, max);
RelRMSE = cell(r, L_max, max);

x_a = cell(r, n);
% x_a_en = cell(r, n);

t_inst = 1:K;
K = t_inst(end);

da_method =@(x_0_en, L, sel) da_seq_bundleIterativeEnsembleKalmanSmoother(x_0_en, y, M, H, R, bIEnKSOptions('OptAlg', @(f, x) QuasiNewtonBFGS(f, x, 1e-6, 50), 'S', 1, 'L', L, 'epsilon', 0.1, 'Localization', sel));;

selector = cell(r, n);

for h=loc
    parfor i1=1:n
        selector{h, i1} = roundSelector(i1, n, h);
    end
end
sel = @(h, i) selector{h, i};

clear
load tmpLOC.mat

tic
for h = loc(1)
    for L = 10 %LWindow
        parfor i = ens
            [ x_a{h, L, i}, x_a_en ] = da_method(x_0_en{i}, L, @(iii) sel(h, iii));
            [ ARMSE(h, L, i), RMSE{h, L, i} ] = averageRootMeanSquareError(x_a{h, L, i}(:,burn_in+1:K), x(:,burn_in+1:K));
            [ ARelRMSE(h, L, i), RelRMSE{h, L, i} ] = averageRelativeRootMeanSquareError(x_a{h, L, i}(:,burn_in+1:K), x(:,burn_in+1:K));
            i
        end
	save tmpLOC
    end
    (h + 1 - loc(1))/(loc(end) - loc(1) + 1)
end
time = toc

%IEnKSLoc_ARMSE = ARMSE;
%clear ARMSE
%IEnKSLoc_RMSE = RMSE;
%clear RMSE
%IEnKSLoc_ARelRMSE = ARelRMSE;
%clear ARelRMSE
%IEnKSLoc_RelRMSE = RelRMSE;
%clear RelRMSE

%save IEnKSOptLocTest

% Optimal: 
% radius = 2 (small ensemble).
% radius = 3 - 4 (big ensemble).
% radius = 5 (very big ensemble).
