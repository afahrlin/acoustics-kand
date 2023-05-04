% Plotting 3D data as scattered grid

% Run in command window
% simname = frequencyHz_key (string)

% TO PLOT YOUR RESULTS
%   - Data should be located one step outside git folder, in a subfolder
%   called /Testdata/simname

% ====================================================

function plot_average(simname)
    
    % Load general data
    location = append('../Testdata/', num2str(simname), '/');
    info = append(location, 'INFO.mat');
    load(info, 'key', 'X_vec', 'Y_vec', 'Z_vec', 'h_t', 'm_t', 'L_x', 'L_y', 'L_z');
    disp('Load done');

    plot_average = true;
    save_average = false;
    load_average = true;
    
    % Plotting parameters
    n = 2;      % Plot every nth point in every direction
    C = 1;      % Size of colored points in scatterplots
    
    % Reshape and select every nth index to visualize
    X = X_vec(1:n:end,1:n:end,1:n:end);
    X = reshape(X, numel(X), 1);
    Y = Y_vec(1:n:end,1:n:end,1:n:end);
    Y = reshape(Y, numel(Y), 1);
    Z = Z_vec(1:n:end,1:n:end,1:n:end);
    Z = reshape(Z, numel(Z), 1);

    

    if save_average
        Average = zeros(size(permute(X_vec, [2,1,3])));
        
        disp(['Number of steps: ', num2str(m_t)])
        for time_step = 1:m_t
            step_file = append(location, num2str(key), '_', num2str(time_step), '.mat');
            load(step_file, 'p');
            
            Average = Average + abs(p);
    
            % Alert every 100th timestep
            if mod(time_step, 100) == 0
                disp([num2str(time_step), '/', num2str(m_t)])
            end
        end
    
        Average = Average/m_t;

        filename = append(location, '/', num2str(key), '_average', '.mat');
        save(filename, 'Average');
    end

    if load_average
        load(append(location, num2str(key), '_average'));
    end

    if plot_average

        Average = permute(Average, [2,1,3]);
        Average = Average(1:n:end, 1:n:end, 1:n:end);
        Average = reshape(Average, numel(Average), 1);

        % Initialize plot ===============================
        % Visualize 4 different angles used a tiledlayout
        t = tiledlayout(7,6,"TileSpacing","compact"); 
        fig = gcf;
        fig.Position = [0, 0, 1000, 1000];
        t.Padding = 'compact';
        title(t,['Average Sound Pressure at Frequency ', simname(1:end-7), ' Hz']);
        cax = [0, 0.2];
    
        % Plot as a scatter plot, Angle 1
        nexttile([3 3]);
        sc1 = scatter3(X, Y, Z, C, Average, 'filled');
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
        sc2 = scatter3(X, Y, Z, C, Average, 'filled');
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
        sc3 = scatter3(X, Y, Z, C, Average, 'filled');
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
        sc4 = scatter3(X, Y, Z, C, Average, 'filled');
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
        drawnow;
        
        % One figure
        figure;
        scatter3(X, Y, Z, C, Average, 'filled');
        title(['Average Sound Pressure at Frequency ', simname(1:end-7), ' Hz']);
        view(-31,74)
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        xlim([0 L_x])
        ylim([0 L_y])
        zlim([0 L_z])
        caxis(cax);
        pbaspect([L_x L_y L_z]);
        cb = colorbar;
        cb.Label.String = 'Sound Pressure';
        drawnow;
    end
end







