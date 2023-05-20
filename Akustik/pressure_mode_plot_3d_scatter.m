% Code to generate plots of mode (n_x,n_y)

L_x = 6.24;
L_y = 3.84;
L_z = 2.4;
n_x = 5;
n_y = 0;
n_z = 0;

f = 343/2*sqrt((n_x/L_x)^2+(n_y/L_y)^2+(n_z/L_z)^2);
disp(['Frequency: ', num2str(f)])

h = 0.077;

m_x = round(L_x/h);
m_y = round(L_y/h);
m_z = round(L_z/h);

x = linspace(0, L_x, m_x);
y = linspace(0, L_y, m_y);
z = linspace(0, L_z, m_z);

[X, Y, Z] = meshgrid(x,y,z);
[X_plot_top, Y_plot_top] = meshgrid(x,y);
[X_plot_side, Z_plot_side] = meshgrid(x,z);
[Y_plot_front, Z_plot_front] = meshgrid(y,z);

c = [-1 1];

P = pressure_mode(1, X, Y, Z, n_x, n_y, n_z, L_x, L_y, L_z);
% P = P + pressure_mode(0.999999999, X, Y, Z, 1, 3, 0, L_x, L_y, L_z)/3;
% P = P + pressure_mode(0.999999999, X, Y, Z, 0, 0, 2, L_x, L_y, L_z)/4;
P = P./max(P,[],'all');

% Plotting parameters n = 1, s = 0.7 for h = 0.077
n = 1;      % Plot every nth point in every direction
s = 0.7;      % Size of colored points in scatterplots

% P = abs(P);

% Reshape and select every nth index to visualize
X = X(1:n:end,1:n:end,1:n:end);
X = reshape(X, numel(X), 1);
Y = Y(1:n:end,1:n:end,1:n:end);
Y = reshape(Y, numel(Y), 1);
Z = Z(1:n:end,1:n:end,1:n:end);
Z = reshape(Z, numel(Z), 1);

disp(['Size P: ', num2str(size(P))])

Average = P(1:n:end, 1:n:end, 1:n:end);
Average = reshape(Average, numel(Average), 1);

% f1 = figure('Position', [1,1,645.6,479.2]);
f1 = figure;
scatter3(X, Y, Z, s, Average, 'filled');
title(['Calculated Pressure Distribution of Mode (', ...
    num2str(n_x), ',', num2str(n_y), ',', num2str(n_z), ')']);
view(-32.226481977574039,38.105776182208906)
% position = [963.4000000000001,242.6,645.5999999999999,479.2];
xlabel('X')
ylabel('Y')
zlabel('Z')
xlim([0 L_x])
ylim([0 L_y])
zlim([0 L_z])
caxis(c);
pbaspect([L_x L_y L_z]);
cb = colorbar;
cb.Label.String = 'Normalized Absolute Sound Pressure';

P = permute(P, [2 1 3]);

f2 = figure;
v = -1:0.01:1;
v2 = [-1,-.75,-.5,-.25,0,.25,.50,.75,1];
[C,h] = contourf(X_plot_top, Y_plot_top, transpose(P(:,:,1)), v);
set(h, 'edgecolor','none');
hold on
[C,h] = contour(X_plot_top, Y_plot_top, transpose(P(:,:,1)), v2, 'k');
clabel(C,h,v2,'FontSize',8,'labelspacing', 1000)
% colormap jet
caxis(c);
%cb1 = colorbar;
%cb1.Label.String = 'Normalized Sound Pressure';
title(['Top View of Calculated Pressure Distribution of Mode (', ...
    num2str(n_x), ',', num2str(n_y), ',', num2str(n_z), ')']);
xlabel('x');
ylabel('y');
pbaspect([L_x L_y, 1]);

f3 = figure;
[C,h] = contourf(X_plot_side, Z_plot_side, transpose(reshape(P(:,1,:), numel(x), numel(z))), v);
%colormap jet
set(h, 'edgecolor','none');
hold on
[C,h] = contour(X_plot_side, Z_plot_side, transpose(reshape(P(:,1,:), numel(x), numel(z))), v2, 'k');
clabel(C,h,v2,'FontSize',8,'labelspacing', 1000)
caxis(c);
% cb2 = colorbar;
% cb2.Label.String = 'Normalized Absolute Sound Pressure';
title(['Side View of Calculated Pressure Distribution of Mode (', ...
    num2str(n_x), ',', num2str(n_y), ',', num2str(n_z), ')']);
xlabel('x');
ylabel('z');
pbaspect([L_x L_z, 1]);

% exportgraphics(f1, 'f1.pdf')
% exportgraphics(f2, 'f2.pdf')
% exportgraphics(f3, 'f3.pdf')

% Initialize plot ===============================
% % Visualize 4 different angles used a tiledlayout
% t = tiledlayout(7,6,"TileSpacing","compact"); 
% fig = gcf;
% fig.Position = [0, 0, 1000, 1000];
% t.Padding = 'compact';
% title(t,['Simulated Pressure Distribution of Mode (', ...
%     num2str(n_x), ',', num2str(n_y), ',', num2str(n_z), ')']);
% 
% % Plot as a scatter plot, Angle 1
% nexttile([3 3]);
% sc1 = scatter3(X, Y, Z, s, Average, 'filled');
% view(-31,14)
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% xlim([0 L_x])
% ylim([0 L_y])
% zlim([0 L_z])
% caxis(c);
% pbaspect([L_x L_y L_z]);
% 
% % Plot as a scatter plot, Angle 2
% nexttile([3 3]);
% sc2 = scatter3(X, Y, Z, s, Average, 'filled');
% view(31,14)
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% xlim([0 L_x])
% ylim([0 L_y])
% zlim([0 L_z])
% caxis(c);
% pbaspect([L_x L_y L_z]);
% 
% % Plot as a scatter plot, Angle 3
% nexttile([4 3]);
% sc3 = scatter3(X, Y, Z, s, Average, 'filled');
% view(-31,74)
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% xlim([0 L_x])
% ylim([0 L_y])
% zlim([0 L_z])
% caxis(c);
% pbaspect([L_x L_y L_z]);
% 
% % Plot as a scatter plot, Angle 4
% nexttile([4 3]);
% sc4 = scatter3(X, Y, Z, s, Average, 'filled');
% view(31,74)
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% xlim([0 L_x])
% ylim([0 L_y])
% zlim([0 L_z])
% caxis(c);
% pbaspect([L_x L_y L_z]);
% 
% % Add colorbar
% cb = colorbar;
% cb.Layout.Tile = 'south';
% cb.Label.String = 'Sound Pressure';
% drawnow;

% % TOP VIEW
% figure;
% scatter3(X, Y, Z, s, Average, 'filled');
% title(['Simulated Pressure Distribution of Mode (', ...
%     num2str(n_x), ',', num2str(n_y), ',', num2str(n_z), ')']);
% % view(-31,74)
% % view(-11.189001050199799,75.442523397394581)
% view(-9.207696632780614,72.04258457960853)
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% xlim([0 L_x])
% ylim([0 L_y])
% zlim([0 L_z])
% caxis(c);
% pbaspect([L_x L_y L_z]);
% 
% 
% 
% % SIDE VIEW
% figure;
% scatter3(X, Y, Z, s, Average, 'filled');
% title(['Simulated Pressure Distribution of Mode (', ...
%     num2str(n_x), ',', num2str(n_y), ',', num2str(n_z), ')']);
% %view(-1.7,3)
% view(-1.983043488570098,6.269456098451739)
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% xlim([0 L_x])
% ylim([0 L_y])
% zlim([0 L_z])
% caxis(c);
% pbaspect([L_x L_y L_z]);
% cb = colorbar;
% cb.Label.String = 'Normalized Absolute Sound Pressure';
% drawnow;

% % Initialize plot ===============================
% % % Visualize 4 different angles used a tiledlayout
% % t = tiledlayout(2,4,"TileSpacing","compact"); 
% % fig = gcf;
% % fig.Position = [0, 0, 1000, 1000];
% % t.Padding = 'compact';
% % title(t,['Simulated Pressure Distribution of Mode (', ...
% %     num2str(n_x), ',', num2str(n_y), ',', num2str(n_z), ')']);
% % 
% % Plot as a scatter plot, Angle 1
% % nexttile([2 2]);
% % sc1 = scatter3(X, Y, Z, s, Average, 'filled');
% % view(-9.207696632780614,72.04258457960853)
% % xlabel('X')
% % ylabel('Y')
% % zlabel('Z')
% % xlim([0 L_x])
% % ylim([0 L_y])
% % zlim([0 L_z])
% % caxis(c);
% % pbaspect([L_x L_y L_z]);
% % 
% % Plot as a scatter plot, Angle 1
% % nexttile([2 2]);
% % sc1 = scatter3(X, Y, Z, s, Average, 'filled');
% % view(-1.983043488570098,6.269456098451739)
% % xlabel('X')
% % ylabel('Y')
% % zlabel('Z')
% % xlim([0 L_x])
% % ylim([0 L_y])
% % zlim([0 L_z])
% % caxis(c);
% % pbaspect([L_x L_y L_z]);



function p = pressure_mode(C, x, y, z, n_x, n_y, n_z, L_x, L_y, L_z)
    p = C*cos(n_x*pi*x/L_x).*cos(n_y*pi*y/L_y).*cos(n_z*pi*z/L_z);
end

