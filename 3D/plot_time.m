% plot single timestep

% Plotting parameters
n = 1;      % Plot every nth point in every direction
C = 1;      % Size of colored points in scatterplots

m_x = 107;
m_y = 65;
m_z = 41;
m = m_x*m_y*m_z;

% Define boundaries (m)
x_l = 0;           % Left boundary of x
x_r = 6.2;            % Right boundary of x
L_x = x_r-x_l;      % Length of x interval
y_l = 0;           % Left boundary of y
y_r = 3.8;            % Right boundary of y
L_y = y_r-y_l;      % Length of y interval
z_l = 0;           % Left boundary of z
z_r = 2.4;            % Right boundary of z
L_z = z_r-z_l;      % Length of z interval

h_x = L_x / (m_x - 1);
x_vec = linspace(x_l, x_r, m_x);
h_y = L_y / (m_y - 1);
y_vec = linspace(y_l, y_r, m_y);
h_z = L_z / (m_z - 1);
z_vec = linspace(z_l, z_r, m_z);
[X_vec, Y_vec, Z_vec] = meshgrid(x_vec, y_vec, z_vec);

X = X_vec(1:n:end,1:n:end,1:n:end);
X = reshape(X, numel(X), 1);
Y = Y_vec(1:n:end,1:n:end,1:n:end);
Y = reshape(Y, numel(Y), 1);
Z = Z_vec(1:n:end,1:n:end,1:n:end);
Z = reshape(Z, numel(Z), 1);

U = zeros(numel(Z),1);

% Initialize plot ===============================
% Visualize 4 different angles used a tiledlayout
t = tiledlayout(7,6,"TileSpacing","compact"); 
fig = gcf;
fig.Position = [0, 0, 1000, 1000];
t.Padding = 'compact';
%title(t,append('Sound Pressure at Time: 0 s. Frequency: ', num2str(f), ' Hz.'));
cax = [-2, 2];

% Plot as a scatter plot, Angle 1
nexttile([3 3]);
sc1 = scatter3(X, Y, Z, C, U, 'filled');
view(-31,14)
xlabel('X')
ylabel('Y')
zlabel('Z')
xlim([0 L_x])
ylim([0 L_y])
zlim([0 L_z])
caxis(cax);
pbaspect([L_x L_y L_z]);

% Plot as a scatter plot, Angle 2
nexttile([3 3]);
sc2 = scatter3(X, Y, Z, C, U, 'filled');
view(31,14)
xlabel('X')
ylabel('Y')
zlabel('Z')
xlim([0 L_x])
ylim([0 L_y])
zlim([0 L_z])
caxis(cax);
pbaspect([L_x L_y L_z]);

% Plot as a scatter plot, Angle 3
nexttile([4 3]);
sc3 = scatter3(X, Y, Z, C, U, 'filled');
view(-31,74)
xlabel('X')
ylabel('Y')
zlabel('Z')
xlim([0 L_x])
ylim([0 L_y])
zlim([0 L_z])
caxis(cax);
pbaspect([L_x L_y L_z]);

% Plot as a scatter plot, Angle 4
nexttile([4 3]);
sc4 = scatter3(X, Y, Z, C, U, 'filled');
view(31,74)
xlabel('X')
ylabel('Y')
zlabel('Z')
xlim([0 L_x])
ylim([0 L_y])
zlim([0 L_z])
caxis(cax);
pbaspect([L_x L_y L_z]);

% Add colorbar
cb = colorbar;
cb.Layout.Tile = 'south';
cb.Label.String = 'Sound Pressure';
% pause(1);

load('/Users/alvafahrlin/Documents/VSCode/Kand/acoustics-kand/3D/uref.mat', 'uref');
size(uref)
U = permute(uref, [2,1,3]);
size(U)

% U = U(1:n:end, 1:n:end, 1:n:end);
U = reshape(U, numel(U), 1);
size(U)

sc1.CData = U;
sc2.CData = U;
sc3.CData = U;
sc4.CData = U;
%title(t,append('Sound Pressure at Time: ', num2str((time_step-1)*h_t), ' s. Frequency: ', num2str(f), ' Hz.'));
drawnow;

