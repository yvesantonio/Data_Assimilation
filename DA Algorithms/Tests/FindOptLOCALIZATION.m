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
%%%%%%%%%%%%%% Run 1:
% % % % % % % min = 5;
% % % % % % % max = 20;
% % % % % % % loc = 1:6;
%%%%%%%%%%%%%% FILE GENERATED: ETKFOptLocTest0
%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Run 2:
min = 21;
max = 35;
loc = 3:5;
%%%%%%%%%%%%%% FILE GENERATED: ETKFOptLocTest1
%%%%%%%%%%%%%%

ens = min:max;
r = length(loc);

ARMSE = zeros(r, max);
RMSE = cell(r, max);
ARelRMSE = zeros(r, max);
RelRMSE = cell(r, max);

x_a = cell(r, n);
% x_a_en = cell(r, n);

t_inst = 1:K;
K = t_inst(end);

da_method = @(x_0_en, selector) da_seq_EnsembleTransformKalmanFilter(x_0_en, y, M, H, R, t_inst, 'Localization', selector);

selector = cell(r, n);

for h=loc
    parfor i1=1:n
        selector{h, i1} = roundSelector(i1, n, h);
    end
end
sel = @(h, i) selector{h, i};

tic
for h = loc
    parfor i=ens
        [ x_a{h, i}, x_a_en ] = da_method(x_0_en{i}, @(iii) sel(h, iii));
        [ ARMSE(h, i), RMSE{h, i} ] = averageRootMeanSquareError(x_a{h, i}(:,burn_in+1:K), x(:,burn_in+1:K));
        [ ARelRMSE(h, i), RelRMSE{h, i} ] = averageRelativeRootMeanSquareError(x_a{h, i}(:,burn_in+1:K), x(:,burn_in+1:K));
        i
    end
    save ETKFOptLocTest1
    (h + 1 - loc(1))/(loc(end) - loc(1) + 1)
end
time = toc

ETKFLoc_ARMSE = ARMSE;
clear ARMSE
ETKFLoc_RMSE = RMSE;
clear RMSE
ETKFLoc_ARelRMSE = ARelRMSE;
clear ARelRMSE
ETKFLoc_RelRMSE = RelRMSE;
clear RelRMSE

save ETKFOptLocTest1

% Optimal: 
% radius = 2 (small ensemble).
% radius = 3 - 4 (big ensemble).
% radius = 5 (very big ensemble).
