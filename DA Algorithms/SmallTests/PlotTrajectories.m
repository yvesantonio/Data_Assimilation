clear
close all
load ETKF

figure
time = 1:0.6*samples_circle;
plot(x(1,time), x(2,time), 'g')

ens = 3;
for i=1:ens
    hold on
    x_plot = squeeze(x_a_en{ens}(1,i, time));
    y_plot = squeeze(x_a_en{ens}(2,i, time));
    plot(x_plot, y_plot, 'b*')
end

figure
time = samples_circle+1:1.9*samples_circle;
plot(x(1,time), x(2,time), 'g')

ens = 3;
for i=1:ens
    hold on
    x_plot = squeeze(x_a_en{ens}(1,i, time));
    y_plot = squeeze(x_a_en{ens}(2,i, time));
    plot(x_plot, y_plot, 'b*')
end


figure
time = 7*samples_circle:8*samples_circle+1;
plot(x(1,time), x(2,time), 'g')

ens = 3;
for i=1:ens
    hold on
    x_plot = squeeze(x_a_en{ens}(1,i, time));
    y_plot = squeeze(x_a_en{ens}(2,i, time));
    plot(x_plot, y_plot, 'b*')
end

figure
time = 20*samples_circle:21*samples_circle+1;
plot(x(1,time), x(2,time), 'g')

ens = 3;
for i=1:ens
    hold on
    x_plot = squeeze(x_a_en{ens}(1,i, time));
    y_plot = squeeze(x_a_en{ens}(2,i, time));
    plot(x_plot, y_plot, 'b*')
end