% ====================================================
% Test of visualization of some 3D data with motion
% ====================================================


function plot_3D()
    % Load data from previous calculations
    filename = '3D_first.mat';
    load(filename, 'X_vec', 'Y_vec', 'Z_vec', 'U', 'h_t', 'm_t', 'L_x', 'L_y', 'L_z');
    disp('Load done');
    U = permute(U,[2,1,3,4]);
    disp('Perm done');
    
    % Initialize video
    Video = VideoWriter('3D_first', 'MPEG-4');
    Video.FrameRate = 60;
    open(Video)
    
    % Plotting parameters
    n = 3;      % Plot every nth point in every direction
    m = 5;      % Plot every mth time step
    C = 1;      % Size of colored points in scatterplots
    
    % Reshape and select every nth index to visualize
    X = X_vec(1:n:end,1:n:end,1:n:end);
    X = reshape(X, numel(X), 1);
    Y = Y_vec(1:n:end,1:n:end,1:n:end);
    Y = reshape(Y, numel(Y), 1);
    Z = Z_vec(1:n:end,1:n:end,1:n:end);
    Z = reshape(Z, numel(Z), 1);
    
    u = U;
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
    sc1 = scatter3(X, Y, Z, C, U, 'filled');
    view(-31,14)
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    pbaspect([L_x L_y L_z]);

    % Plot as a scatter plot, Angle 2
    nexttile([3 3]);
    sc2 = scatter3(X, Y, Z, C, U, 'filled');
    view(31,14)
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    pbaspect([L_x L_y L_z]);

    % Plot as a scatter plot, Angle 3
    nexttile([4 3]);
    sc3 = scatter3(X, Y, Z, C, U, 'filled');
    view(-31,74)
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    %ylim([y_l y_r]);
    pbaspect([L_x L_y L_z]);

    % Plot as a scatter plot, Angle 4
    nexttile([4 3]);
    sc4 = scatter3(X, Y, Z, C, U, 'filled');
    view(31,74)
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    %ylim([y_l y_r]);
    pbaspect([L_x L_y L_z]);
    
    % Add colorbar
    cb = colorbar;
    caxis([-0.1,0.1]);
    cb.Layout.Tile = 'south';
    cb.Label.String = 'Sound Pressure';
    pause(1);
    
    for i = 1:m_t/m
        % Reshape next time step
        U = u(:,:,:,i*m);
        U = U(1:n:end,1:n:end,1:n:end);
        U = reshape(U, numel(U), 1);
        
        % Update data and draw it by updating title
        sc1.CData = U;
        sc2.CData = U;
        sc3.CData = U;
        sc4.CData = U;
        title(t,['Visualization! Time: ', num2str(h_t*(i*m-1))]);
        
        % Write the timestep as a frame to the video
        frame = getframe(gcf);
        writeVideo(Video, frame);
    end
    
    % Close video
    close(Video)
end







