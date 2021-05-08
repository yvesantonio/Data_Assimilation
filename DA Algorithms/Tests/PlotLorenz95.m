% Plot Lorenz Realizations
clear
close all
load BookTestData.mat

x = x(:, 2000:2999);
K = 1000;

x_dis = [ x; x(1,:) ]; % Circular System: fancy looking plot.
[sates_grid, t_grid] = meshgrid(0:n, 1:K);

figure('Renderer', 'painters', 'Position', [300 200 560 317])
surface(t_grid, sates_grid, x_dis')
colormap(fliplr(gray()')')
title('States evolution')
xlabel('time');
ylabel('State Index');
zlabel('State');
shading interp
colorbar
view(30, 60)

figure('Renderer', 'painters', 'Position', [300 200 560 317])
surface(t_grid, sates_grid, x_dis')
colormap(fliplr(gray()')')
title('States evolution')
xlabel('time');
ylabel('State Index');
zlabel('State');
shading interp
colorbar
view(0, 90)