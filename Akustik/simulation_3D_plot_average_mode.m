% 3D simulation of the acoustic wave equation
% 
% FINAL CODE
% 
% Simulation of the acoustic wave equation in
% three dimensions. Plots sound pressure over 
% time. 
% 
% TO SAVE YOUR RESULTS
%   - Create a folder one step outside acoustics-kand named Testdata
%   - Your tests will be saved there (if true) in a subfolder called
%     frequencyHz_key, where key is a random 4 digit nr
% 
% =========================================================================

function simulation_3D_plot_average_mode()

    % =====================================================================
    % Model parameters
    f = 137;
    w = 2*pi*f;     % Angular frequency
    n_x = 4;
    n_y = 0;
    n_z = 1;
    
    T = 1/f*20+1/f/2;
    
    % Define boundaries (m)
    x_l = 0;            % Left boundary of x
    x_r = 6.24;         % Right boundary of x
    L_x = x_r-x_l;      % Length of x interval
    y_l = 0;            % Left boundary of y
    y_r = 3.84;         % Right boundary of y
    L_y = y_r-y_l;      % Length of y interval
    z_l = 0;            % Left boundary of z
    z_r = 2.4;          % Right boundary of z
    L_z = z_r-z_l;      % Length of z interval
    
    % Number of grid points
%     m_x = 185;      % High-res grid
%     m_y = 113;
%     m_z = 71;
    m_x = 81;       % Low-res grid
    m_y = 49;
    m_z = 31;
%     m_x = 29;       % Very low-res grid
%     m_y = 17;
%     m_z = 11;
%     m_x = 13;       % Eigenvalue tests
%     m_y = 9;
%     m_z = 5;
%     m_x = 19;
%     m_y = 11;
%     m_z = 7;
%     m_x = 625;       % UPPMAX test 1
%     m_y = 385;
%     m_z = 241;
    m = m_x*m_y*m_z;
    
    % =====================================================================
    % PDE parameters
    
    c = 343;            % Wave speed (m/s)
    beta_2 = c;
    
    % Determine beta_3 based on reflection coefficient
    R = 1;              % Reflection rate, must be between 0 and 1
    p = [-0.454, 1.34, -1.883, 1-R];
    r = roots(p);
    beta_3 = r(3);
    if R == 1
        beta_3 = 0;
    end
    
    % =====================================================================
    % Discretization
    
    % Spatial discretization
    h_x = L_x / (m_x-1);
    x_vec = linspace(x_l, x_r, m_x);
    h_y = L_y / (m_y-1);
    y_vec = linspace(y_l, y_r, m_y);
    h_z = L_z / (m_z-1);
    z_vec = linspace(z_l, z_r, m_z);
    [X_vec, Y_vec, Z_vec] = meshgrid(x_vec, y_vec, z_vec);
    
    disp(h_x);
    disp(h_y);
    disp(h_z);
    
    % Time discretization
%     h_t = 0.1*max([h_x, h_y, h_z])/c;
    h = max([h_x, h_y, h_z]);
    if R < 0.2
%         h_t = 0.001*max([h_x, h_y, h_z]);
        h_t = 0.95*2.6e-4*(1.556*h^2+3.866*h);
    else
%         h_t = 0.002*max([h_x, h_y, h_z]);
        h_t = 0.95*1.85*2.6e-4*(1.556*h^2+3.866*h);
    end
    m_t = round(T/h_t,0);
    h_t = T/m_t;
    
    disp('Discretization Done')
    % =====================================================================
    % Point source parameters
    
    % Point source amplitude
    vol = 80;                           % Speaker volume
    Pvol = 20*10^(-6) * 10^(vol/20);    % Pressure 1m from speaker
    K = 17.0852;
    amp = Pvol/(K*h_y*h_z);        % Amplitude of point source
    amp_ps = amp;                       % Constant amp of point source

    % ================================================================
    % SBP-SAT method
    
    % Get D2 operator - x
    [~, HI_x, ~, D2_x, e_lx, e_rx, d1_lx, d1_rx] = sbp_cent_4th(m_x, h_x);
    % SBP-SAT
    D_x = c^2*(D2_x + HI_x*e_lx'*d1_lx - HI_x*e_rx'*d1_rx);
    E_x = -c^2*beta_3/beta_2*HI_x*(e_lx'*e_lx + e_rx'*e_rx);

    % Get D2 operator - y
    [~, HI_y, ~, D2_y, e_ly, e_ry, d1_ly, d1_ry] = sbp_cent_4th(m_y, h_y);
    % SBP-SAT
    D_y = c^2*(D2_y + HI_y*e_ly'*d1_ly - HI_y*e_ry'*d1_ry);
    E_y = -c^2*beta_3/beta_2*HI_y*(e_ly'*e_ly + e_ry'*e_ry);
    
    % Get D2 operator - z
    [~, HI_z, ~, D2_z, e_lz, e_rz, d1_lz, d1_rz] = sbp_cent_4th(m_z, h_z);
    % SBP-SAT
    D_z = c^2*(D2_z + HI_z*e_lz'*d1_lz - HI_z*e_rz'*d1_rz);
    E_z = -c^2*beta_3/beta_2*HI_z*(e_lz'*e_lz + e_rz'*e_rz);
    
    % SBP operators
    D_xy = sparse(kron(speye(m_y), D_x) + kron(D_y, speye(m_x)));
    D = sparse(kron(speye(m_z), D_xy) + kron(D_z, speye(m_x*m_y)));
    disp('D-Operator Done')
    
    E_xy = sparse(kron(speye(m_y), E_x) + kron(E_y, speye(m_x)));
    E = sparse(kron(speye(m_z), E_xy) + kron(E_z, speye(m_x*m_y)));
    disp('E-Operator Done')
    
    % Construct matrix A: v_t = Av with v = [p, p_t]^T
    % [0, I;
    %  D, E]
    A = [sparse(m,m), speye(m); D, E];
    disp('A-Matrix Done')
    
    % Set initial values
    u = zeros(2*m, 1);      % Initial pressure deviation
    Average = zeros(m,1);
    t = 0;                  % Initial time
    
    % =====================================================================
    % INFOSTRING
    disp(['Frequency: ', num2str(w/(2*pi)), ' Hz']);
    disp(['Number of gridpoints: ', num2str(size(X_vec))])
    disp(['Simulation time: ', num2str(T), 's'])
    disp(['Number of steps: ', num2str(m_t)])
    
    % =====================================================================
    % Plot and time step
    % Step through time with RK4
    for time_step = 1:m_t
        [u,t] = steprk4(u, t, h_t);
        Average = Average + abs(u(1:m));
        
        % Alert every 100th timestep
        if mod(time_step, 100) == 0
            disp(time_step)
        end
    end
    
%     Average = Average./m_t;

    P = u(1:m);
    P = reshape(P, m_x, m_y, m_z);
    
    c = [-1 1];
    % Plotting parameters n = 1, s = 0.7 for h = 0.077
    n = 1;      % Plot every nth point in every direction
    s = 0.7;      % Size of colored points in scatterplots
    
    % Reshape and select every nth index to visualize
    X = X_vec(1:n:end,1:n:end,1:n:end);
    X = reshape(X, numel(X), 1);
    Y = Y_vec(1:n:end,1:n:end,1:n:end);
    Y = reshape(Y, numel(Y), 1);
    Z = Z_vec(1:n:end,1:n:end,1:n:end);
    Z = reshape(Z, numel(Z), 1);
    
    [X_plot_top, Y_plot_top] = meshgrid(x_vec,y_vec);
    [X_plot_side, Z_plot_side] = meshgrid(x_vec,z_vec);
    P = permute(P, [2 1 3]);
    if abs(min(P,[],'all')) > abs(max(P,[],'all'))
        P = P./abs(min(P,[],'all'));
    else
        P = P./abs(max(P,[],'all'));
    end
    disp(['Max value: ', num2str(max(P,[],'all'))])
    disp(['Min value: ', num2str(min(P,[],'all'))])
    
    disp(['Size P: ', num2str(size(P))])
    
    Average = P(1:n:end, 1:n:end, 1:n:end);
    Average = reshape(Average, numel(Average), 1);

    % f1 = figure('Position', [1,1,645.6,479.2]);
    f1 = figure;
    scatter3(X, Y, Z, s, Average, 'filled');
    title(['Simulated Pressure Distribution of Mode (', ...
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
    title(['Top View of Simulated Pressure Distribution of Mode (', ...
        num2str(n_x), ',', num2str(n_y), ',', num2str(n_z), ')']);
    xlabel('x');
    ylabel('y');
    pbaspect([L_x L_y, 1]);

    f3 = figure;
    [C,h] = contourf(X_plot_side, Z_plot_side, transpose(reshape(P(:,1,:), numel(x_vec), numel(z_vec))), v);
    %colormap jet
    set(h, 'edgecolor','none');
    hold on
    [C,h] = contour(X_plot_side, Z_plot_side, transpose(reshape(P(:,1,:), numel(x_vec), numel(z_vec))), v2, 'k');
    clabel(C,h,v2,'FontSize',8,'labelspacing', 1000)
    caxis(c);
    % cb2 = colorbar;
    % cb2.Label.String = 'Normalized Absolute Sound Pressure';
    title(['Side View of Simulated Pressure Distribution of Mode (', ...
        num2str(n_x), ',', num2str(n_y), ',', num2str(n_z), ')']);
    xlabel('x');
    ylabel('z');
    pbaspect([L_x L_z, 1]);

    % exportgraphics(f1, 'f1.pdf')
    % exportgraphics(f2, 'f2.pdf')
    % exportgraphics(f3, 'f3.pdf')
    
    % =====================================================================
    % Define functions used in code 
    
    % Define rhs of the semi-discrete approximation
    function u_t = rhs(t, u) 
        u_t = A*F2(t, u);
    end
    
    % Update value of point sources
    function v = F2(t, v)
        % For a 'smooth' start
        if t < 0.01
            amp_ps = amp*t/0.01;
        else
            amp_ps = amp;
        end
        
        % One point source in middle of left x-boundary
%         v(round(m_x*m_y*m_z/2, 0)-m_x*round(0.5*m_y, 0)+...
%             round(0.5*m_x*m_y)) = amp_ps*sin(w*t);
        
        % Two point sources on left x-boundary
        v((round(m/2, 0)-m_x*round(0.25*m_y, 0))+round(0.5*m_x*m_y)) = ...
            amp_ps*sin(w*t);
        v((round(m/2, 0)-m_x*round(0.75*m_y, 0))+round(0.5*m_x*m_y)) = ...
            amp_ps*sin(w*t);
    end

    % Time step with rk4
    function [v, t] = steprk4(v, t, dt)
        % Rk4 coefficients
        k1 = dt*rhs(t, v);
        k2 = dt*rhs(t+0.5*dt, v+0.5*k1);
        k3 = dt*rhs(t+0.5*dt, v+0.5*k2);
        k4 = dt*rhs(t+dt, v+k3);

        % Calc new value and time
        v = v + 1/6*(k1 + 2*k2 + 2*k3 + k4);
        t = t + dt;
        
        % Update point sources for the plot
        v = F2(t, v);
    end
    
    % Get reference line of sound propagation for 1D plot
    function a = amplitude(r, Pvol, R)
        a = (1+R)*Pvol./r;
    end
end





