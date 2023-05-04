% Plotting 3D data as scattered grid

% Run in command window
% simname = frequencyHz_key (string)

% TO PLOT YOUR RESULTS
%   - Data should be located one step outside git folder, in a subfolder
%   called /Testdata/simname

% ====================================================

function plot_3D(simname)
    
    % Load general data
    location = append('../Testdata/', num2str(simname), '/');
    info = append(location, 'INFO.mat');
    load(info, 'key', 'X_vec', 'Y_vec', 'Z_vec', 'h_t', 'm_t', 'L_x', 'L_y', 'L_z', 'f');
    disp('Load done');
    
    % Initialize video
    Video = VideoWriter(append(location, num2str(simname), '_video'), 'MPEG-4');
    Video.FrameRate = 60;
    open(Video)
    
    % Plotting parameters
    n = 2;      % Plot every nth point in every direction
    m = 5;      % Plot every mth time step
    C = 1;      % Size of colored points in scatterplots

    % Reshape and select every nth index to visualize
    X = X_vec(1:n:end,1:n:end,1:n/2:end);
    X = reshape(X, numel(X), 1);
    Y = Y_vec(1:n:end,1:n:end,1:n/2:end);
    Y = reshape(Y, numel(Y), 1);
    Z = Z_vec(1:n:end,1:n:end,1:n/2:end);
    Z = reshape(Z, numel(Z), 1);
    
    U = zeros(numel(Z),1);
    
    % Initialize plot ===============================
    % Visualize 4 different angles used a tiledlayout
    t = tiledlayout(7,6,"TileSpacing","compact"); 
    fig = gcf;
    fig.Position = [0, 0, 1000, 1000];
    t.Padding = 'compact';
    title(t,append('Sound Pressure at Time: 0 s. Frequency: ', num2str(f), ' Hz.'));
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
    pause(1);
    
    for time_step = 1:m:m_t
        step_file = append(location, num2str(key), '_', num2str(time_step), '.mat');
        load(step_file, 'p');
        U = permute(p, [2,1,3]);
        
        U = U(1:n:end, 1:n:end, 1:n/2:end);
        U = reshape(U, numel(U), 1);

        sc1.CData = U;
        sc2.CData = U;
        sc3.CData = U;
        sc4.CData = U;
        title(t,append('Sound Pressure at Time: ', num2str((time_step-1)*h_t), ' s. Frequency: ', num2str(f), ' Hz.'));
        drawnow;
        
        % Write the timestep as a frame to the video
        frame = getframe(gcf);
        writeVideo(Video, frame);
    end
    
    % Close video
    close(Video)
end







