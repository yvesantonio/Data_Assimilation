% surface plot Inflation
clear
close all

% load ETKFOptInflTest0.mat
% load ETKFOptInflTest1.mat
% load ETKFOptInflTest2.mat
% load ETKFOptInflTest3.mat

[ infl_grid, ens_grid ] = meshgrid(infl, ens);

figure('Renderer', 'painters', 'Position', [500 300 560 317])
surface(infl_grid, ens_grid, ETKFInfl_ARMSE(:,ens)')
title('ARMSE')
xlabel('Inflation');
ylabel('Ensemble Size');
zlabel('ARMSE');
shading interp
%colorbar
view(0,90)

figure('Renderer', 'painters', 'Position', [500 300 560 317])
surface(infl_grid, ens_grid, ETKFInfl_ARMSE(:,ens)')
title('ARMSE')
xlabel('Inflation');
ylabel('Ensemble Size');
zlabel('ARMSE');
shading interp
colorbar
view(10,-80)