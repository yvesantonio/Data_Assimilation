% surface plot Inflation
clear
close all

load tmpLOC.mat
% load ETKFOptLocTest1.mat

loc = [5 7 10 12];
[ loc_grid, ens_grid ] = meshgrid(loc, ens);

L = 1;

figure
surface(loc_grid, ens_grid, squeeze(ARMSE(loc, L, ens))')
title(['IEnKS L = ' num2str(L) ' with Localization, ARMSE'])
xlabel('Localization Radius');
ylabel('Ensemble Size');
zlabel('ARMSE');
shading interp
colorbar
view(0,90)

L = 5;

figure
surface(loc_grid, ens_grid, squeeze(ARMSE(loc, L, ens))')
title(['IEnKS L = ' num2str(L) ' with Localization, ARMSE'])
xlabel('Localization Radius');
ylabel('Ensemble Size');
zlabel('ARMSE');
shading interp
colorbar
view(0,90)

L = 10;

figure
surface(loc_grid, ens_grid, squeeze(ARMSE(loc, L, ens))')
title(['IEnKS L = ' num2str(L) ' with Localization, ARMSE'])
xlabel('Localization Radius');
ylabel('Ensemble Size');
zlabel('ARMSE');
shading interp
colorbar
view(0,90)