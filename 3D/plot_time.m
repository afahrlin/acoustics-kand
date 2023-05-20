% plot single timestep

function plot_time()
    
    % Plotting parameters
    n = 6;      % Plot every nth point in every direction
    C = 1;      % Size of colored points in scatterplots

    % Import data
    load('Konvergens/t6/INFO.mat', 'X_vec', 'Y_vec', 'Z_vec', 'L_x', 'L_y', 'L_z', 'T', 'f');
    disp('Load info done');
    
    load('Konvergens/t6/u.mat', 'u');
    disp('Load uref done')
    
    time = T;

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
    title(t,append('Sound Pressure at Time: ', num2str(time), ' s. Frequency: ', num2str(f), ' Hz.'));
    cax = [0, 0.1];

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
    
    U = permute(u, [2,1,3]);
    U = U(1:n:end, 1:n:end, 1:n:end);
    %U = 20 * log10(abs(reshape(U, numel(U), 1))./(20*10^(-6)));
    U = reshape(U, numel(U), 1);

    sc1.CData = U;
    sc2.CData = U;
    sc3.CData = U;
    sc4.CData = U;
    
    drawnow;
end
