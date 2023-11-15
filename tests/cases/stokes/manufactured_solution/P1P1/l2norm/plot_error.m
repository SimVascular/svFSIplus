clear all; clc; close all;

figure('units','normalized','outerposition',[0.65 0.4 0.3 0.48]);
hold on;

%% Plot L2 error in velocity
fname = 'l2norm_v.dat';
fid = fopen(fname, 'r');
C = textscan(fid, '%f %f %f %f', 'HeaderLines', 1);
fclose(fid);

dx = C{1};
relEr = C{2};
errL2 = C{3};
[p,~] = polyfit(log(dx), log(relEr), 1);
s1 = sprintf('slope (v) = %.2f', p(1));

h(1) = plot(dx, relEr, 'r-s', 'LineWidth', 2,...
    'MarkerSize', 10, 'MarkerEdgeColor', 'r');
h(2) = plot(dx, errL2, 'r--', 'LineWidth', 2);


%% Plot L2 error in pressure
fname = 'l2norm_p.dat';
fid = fopen(fname, 'r');
C = textscan(fid, '%f %f %f %f', 'HeaderLines', 1);
fclose(fid);

dx = C{1};
relEr = C{2};
errL2 = C{3};
[p,~] = polyfit(log(dx), log(relEr), 1);
s2 = sprintf('slope (p) = %.2f', p(1));

h(3) = plot(dx, relEr, 'k-d', 'LineWidth', 2,...
    'MarkerSize', 10, 'MarkerEdgeColor', 'k');
h(4) = plot(dx, errL2, 'k--', 'LineWidth', 2);


%% General plot features
axis([0.002  0.2 1.0e-4 0.5]);
set(gca, 'FontName', 'Times', 'FontSize', 14, 'LineWidth', 1);
hXLabel = xlabel('$\Delta$ x','interpreter','latex');
hYLabel = ylabel('$||u-U_e||_2 / ||U_e||_2$','interpreter','latex');
set([hXLabel, hYLabel], 'FontName', 'Times', 'FontSize', 20,...
    'FontWeight', 'bold');
set(gca, 'XScale', 'log', 'YScale', 'log');
set( gca, 'Box', 'on', 'TickDir'     , 'out', ...
    'TickLength'  , [.01 .01], ...
    'XMinorTick'  , 'on', ...
    'YMinorTick'  , 'on', ...
    'YGrid'       , 'on', ...
    'XGrid'       , 'on', ...
    'XColor'      , [0 0 0 ], ...
    'YColor'      , [0 0 0 ], ...
    'LineWidth'   , 1 );

hleg = legend(h([1 3]), 'v', 'p', 'location', 'NorthWest');
set(hleg, 'FontSize', 15);
leg_pos = get(hleg,'position');
set(hleg, 'position', [leg_pos(1) leg_pos(2) 1.5*leg_pos(3) leg_pos(4)]);
hleg.ItemTokenSize = [50 50];

dim = [0.6 0.2 0.25 0.1];
str = {'P1P1 (stab)', s1, s2};
annotation('textbox',dim,'String',str,'FitBoxToText','on',...
    'FontSize',15, 'BackgroundColor','w');

hold off;

% EOF