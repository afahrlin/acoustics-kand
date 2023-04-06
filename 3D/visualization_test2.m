% ====================================================
% Test of visualization of some 3D data with motion
% ====================================================


function visualization_test2()
    filename = 'test.mat';
    load(filename, 'X_vec', 'Y_vec', 'Z_vec', 'u', 'h_t', 'm_t', 'L_x', 'L_y', 'L_z');
    
    % Reshape and select every nth index to visualize
    n = 3;
    X = X_vec(1:n:end,1:n:end,1:n:end);
    X = reshape(X, numel(X), 1);
    Y = Y_vec(1:n:end,1:n:end,1:n:end);
    Y = reshape(Y, numel(Y), 1);
    Z = Z_vec(1:n:end,1:n:end,1:n:end);
    Z = reshape(Z, numel(Z), 1);
    
    U = u(:,:,:,1);
    U = U(1:n:end,1:n:end,1:n:end);
    U = reshape(U, numel(U), 1);
    
    % Initialize plot ===============================
    % Visualize 4 different angles used a tiledlayout
    t = tiledlayout(7,6,"TileSpacing","compact");
    fig = gcf;
    fig.Position = [0, 0, 1000, 1000];
    t.Padding = 'compact';
    title(t,'Visualization 2! Time: 0');

    % Plot as a scatter plot, Angle 1
    nexttile([3 3]);
    sc1 = scatter3(X, Y, Z, 1, U, 'filled');
    view(-31,14)
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    pbaspect([L_x L_y L_z]);

    % Plot as a scatter plot, Angle 2
    nexttile([3 3]);
    sc2 = scatter3(X, Y, Z, 1, U, 'filled');
    view(31,14)
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    pbaspect([L_x L_y L_z]);

    % Plot as a scatter plot, Angle 3
    nexttile([4 3]);
    sc3 = scatter3(X, Y, Z, 1, U, 'filled');
    view(-31,74)
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    %ylim([y_l y_r]);
    pbaspect([L_x L_y L_z]);

    % Plot as a scatter plot, Angle 4
    nexttile([4 3]);
    sc4 = scatter3(X, Y, Z, 1, U, 'filled');
    view(31,74)
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    %ylim([y_l y_r]);
    pbaspect([L_x L_y L_z]);
    
    % Add colorbar
    cb = colorbar;
    cb.Layout.Tile = 'south';
    cb.Label.String = 'Sound Pressure';
    pause(1);
    
    m = 1;
    
    for i = 1:m_t/m
        tic
        U = u(:,:,:,i*m);
        U = U(1:n:end,1:n:end,1:n:end);
        U = reshape(U, numel(U), 1);
        toc
        tic
        sc1.CData = U;
        sc2.CData = U;
        sc3.CData = U;
        sc4.CData = U;
        toc
        tic
        title(t,['Visualization 2! Time: ', num2str(h_t*(i*m-1))]);
        pause(0.001);
        toc
    end
end







