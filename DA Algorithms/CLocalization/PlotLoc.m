% surface plot Inflation
clear
close all

rho = roundRho(40, 8);
x = 1:40;
[ x_grid y_grid ] = meshgrid(x, x);
figure('Renderer', 'painters', 'Position', [500 300 400*4/3 400 ])
set(gca,'Ydir','reverse')
surface(x_grid, y_grid, rho)
title('Rho, c=8, Gaspari-Cohn function')
xlabel('x');
ylabel('y');
zlabel('rho(x,y)');
shading interp
colormap(jet())
colorbar
view(0,90)


load CLocETKFTest
[ loc_grid, ens_grid ] = meshgrid(loc, ens);

figure('Renderer', 'painters', 'Position', [500 300 560 317])
surface(loc_grid, ens_grid, CLocETKF_ARelRMSE(loc,ens)')
title('ETKF, Relative ARMSE')
xlabel('Localization Radius, c');
ylabel('Ensemble Size');
zlabel('ARMSE');
shading interp
colorbar
view(0,90)

% % % % % % % % % % % % % % % % % % figure('Renderer', 'painters', 'Position', [500 300 560 317])
% % % % % % % % % % % % % % % % % % surface(loc_grid, ens_grid, CLocEnKF_ARelRMSE(loc,ens)')
% % % % % % % % % % % % % % % % % % title('ETKF, Relative ARMSE')
% % % % % % % % % % % % % % % % % % xlabel('Localization Radius, c');
% % % % % % % % % % % % % % % % % % ylabel('Ensemble Size');
% % % % % % % % % % % % % % % % % % zlabel('ARMSE');
% % % % % % % % % % % % % % % % % % shading interp
% % % % % % % % % % % % % % % % % % colorbar
% % % % % % % % % % % % % % % % % % view(0,0)

load CLocEnKFTest

figure('Renderer', 'painters', 'Position', [500 300 560 317])
surface(loc_grid, ens_grid, CLocEnKF_ARelRMSE(loc,ens)')
title('EnKF, Relative ARMSE')
xlabel('Localization Radius, c');
ylabel('Ensemble Size');
zlabel('ARMSE');
shading interp
colorbar
view(0,90)