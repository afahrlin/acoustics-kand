% ====================================================
% Test of visualization of some 3D data
% 
% Some different ways of visualizing 3D data. 
% Pay attention to the titles of the figures to know
% which viewing angles go with which visualization
% method. 
% 
% Alternativ 2 var krånglig att få till med tiles. 
% ====================================================

function visualization_tests()
    % Define boundaries (m)
    x_l = 0;            % Left boundary of x
    x_r = 10;           % Right boundary of x
    L_x = x_r-x_l;      % Length of x interval
    
    y_l = 0;            % Left boundary of y
    y_r = 5;            % Right boundary of y
    L_y = y_r-y_l;      % Length of y interval
    
    z_l = 0;            % Left boundary of z
    z_r = 5;            % Right boundary of z
    L_z = z_r-z_l;      % Length of z interval
    
    % Used in initial data f0
    x_0 = L_x/2;
    y_0 = L_y/2;
    z_0 = L_z/2;
    
    % Number of grid points
    m_x = 200;
    m_y = 100;
    m_z = 100;
    m = m_x*m_y*m_z;
    
    % Spatial discretization
    h_x = L_x / (m_x - 1);
    x_vec = linspace(x_l, x_r, m_x);
    
    h_y = L_y / (m_y - 1);
    y_vec = linspace(y_l, y_r, m_y);
    
    h_z = L_z / (m_z - 1);
    z_vec = linspace(z_l, z_r, m_z);
    
    % Create meshgrid
    [X_vec, Y_vec, Z_vec] = meshgrid(x_vec, y_vec, z_vec);
    
    % Get initial values
    u = f0(X_vec, Y_vec, Z_vec);
    
    % ======================================================
    % Visualisation
    plot_alternative = 2;
    
    % ALTERNATIVE 1 ========================================
    % Reshape and select every nth index to visualize
    if plot_alternative == 1
        n = 27;
        X = reshape(X_vec, m, 1);
        X = X(1:n:end,1);
        Y = reshape(Y_vec, m, 1);
        Y = Y(1:n:end,1);
        Z = reshape(Z_vec, m, 1);
        Z = Z(1:n:end,1);
        U = reshape(u, m, 1);
        U = U(1:n:end,1);

        % Plot as a scatter plot
        figure;
        scatter3(X, Y, Z, 1, U, 'filled')
        view(-31,14)
        title('Visualization 1')
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        pbaspect([L_x L_y L_z]);
        cb = colorbar;
        cb.Label.String = 'Sound Pressure';
    end
    
    % ALTERNATIVE 2 ========================================
    % Reshape and select every nth index to visualize
    if plot_alternative == 2
        n = 3;
        X = X_vec(1:n:end,1:n:end,1:n:end);
        X = reshape(X, numel(X), 1);
        Y = Y_vec(1:n:end,1:n:end,1:n:end);
        Y = reshape(Y, numel(Y), 1);
        Z = Z_vec(1:n:end,1:n:end,1:n:end);
        Z = reshape(Z, numel(Z), 1);
        U = u(1:n:end,1:n:end,1:n:end);
        U = reshape(U, numel(U), 1);

        % Visualize 4 different angles used a tiledlayout
        t = tiledlayout(7,6,"TileSpacing","compact");
        fig = gcf;
        fig.Position = [0, 0, 1000, 1000];
        t.Padding = 'compact';
        title(t,'Visualization 2');
        
        % Plot as a scatter plot, Angle 1
        ax1 = nexttile([3 3]);
        scatter3(X, Y, Z, 1, U, 'filled')
        view(-31,14)
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        pbaspect([L_x L_y L_z]);
        

        % Plot as a scatter plot, Angle 2
        ax2 = nexttile([3 3]);
        scatter3(X, Y, Z, 1, U, 'filled')
        view(31,14)
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        pbaspect([L_x L_y L_z]);

        % Plot as a scatter plot, Angle 3
        ax3 = nexttile([4 3]);
        scatter3(X, Y, Z, 1, U, 'filled')
        view(-31,76)
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        ylim([y_l y_r]);
        pbaspect([L_x L_y L_z]);

        % Plot as a scatter plot, Angle 3
        ax4 = nexttile([4 3]);
        scatter3(X, Y, Z, 1, U, 'filled')
        view(31,76)
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        ylim([y_l y_r]);
        pbaspect([L_x L_y L_z]);

        cb = colorbar;
        cb.Layout.Tile = 'south';
        cb.Label.String = 'Sound Pressure';
    end
    
    % ALTERNATIVE 3 ========================================
    % Reshape and select every nth index to visualize
    if plot_alternative == 3
        n = 18;
        X = X_vec(1:n:end,1:n:end,1:n:end);
        X = reshape(X, numel(X), 1);
        Y = Y_vec(1:n:end,1:n:end,1:n:end);
        Y = reshape(Y, numel(Y), 1);
        Z = Z_vec(1:n:end,1:n:end,1:n:end);
        Z = reshape(Z, numel(Z), 1);
        U = u(1:n:end,1:n:end,1:n:end);
        U = reshape(U, numel(U), 1);

        % Plot as a scatter plot
        figure;
        scatter3(X, Y, Z, 20, U, 'filled')
        view(-31,14)
        title('Visualization 3')
        xlabel('X')
        ylabel('Y')
        zlabel('Z')
        pbaspect([L_x L_y L_z]);
        cb = colorbar;
        cb.Label.String = 'Sound Pressure';
    end
    
    
    % Function to visualize in 3d
    function v = f0(X, Y, Z)
        % Some different functions to test visualize
        % v = 0*X;
        % v = sin(X);
        % v = sin(X+Y+Z);
        % v = sin(X.*Y.*Z);
        
        r = sqrt((X-x_0).^2+(Y-y_0).^2+(Z-z_0).^2);
        v = sin(2*r);
    end
    
end







