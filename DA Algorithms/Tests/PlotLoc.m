% surface plot Inflation
clear
close all

load ETKFOptLocTest0.mat
% load ETKFOptLocTest1.mat
[ loc_grid, ens_grid ] = meshgrid(loc, ens);

figure
surface(loc_grid, ens_grid, ETKFLoc_ARMSE(loc,ens)')
title('ARMSE')
xlabel('Localization Radius');
ylabel('Ensemble Size');
zlabel('ARMSE');
shading interp
colorbar
view(0,90)

figure
surface(loc_grid, ens_grid, ETKFLoc_ARMSE(loc,ens)')
title('ARMSE')
xlabel('Localization Radius');
ylabel('Ensemble Size');
zlabel('ARMSE');
shading interp
colorbar
view(0,0)