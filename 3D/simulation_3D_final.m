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

function simulation_3D_final()
    calc_eigenvalues = false;       % If true, calc eigenvalues of A
    plot_time_steps_1d = false;     % If true, plot time steps 1d
    plot_time_steps_2d = true;     % If true, plot time steps 2d
    save_time_steps = false;        % If true, save time steps
    
    % =====================================================================
    % Model parameters
    
    T = 1;              % Final time (seconds)
    f = 110;            % Frequency (Hz)
    s = 1;              % plot every s time-steps
    
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
    if R == 1           % When the root is 0, it has a different index and
        beta_3 = 0;     % must be set manually. 
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
    
    disp(['Step size x: ', num2str(h_x)]);
    disp(['Step size y: ', num2str(h_y)]);
    disp(['Step size z: ', num2str(h_z)]);
    
    % Time discretization
    h = max([h_x, h_y, h_z]);   % Since h is not exactly the same in xyz
    % Stability constraints derived in report, stable for all tests so far
    if R < 0.2
        h_t = 0.95*2.6e-4*(1.556*h^2+3.866*h);
    else
        h_t = 0.95*1.85*2.6e-4*(1.556*h^2+3.866*h);
    end
    m_t = round(T/h_t,0);
    h_t = T/m_t;
    
    disp('Discretization Done')
    % =====================================================================
    % Point source parameters
    
    % Point source frequency
%     lambda = max([L_x L_y L_z]);    % Longest wave resonant with the room
%     k = 4;                          % Which overtone
%     f = k*c/lambda;                 % Frequency
    w = 2*pi*f;     % Angular frequency
    
    % Point source amplitude
    vol = 80;                       % Speaker volume (dB 1m from speaker)
    Pvol = 20*10^(-6)*10^(vol/20);  % Pressure 1m from speaker
    K = 17.0852;                    % From test with a very fine grid, h=0.01
    amp = Pvol/(K*h_y*h_z);         % Scaling of amplitude in point source
    amp_ps = amp;                   % Constant amp of point source
    
    % Plot amplitude
    z = [-Pvol*4 Pvol*4];

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
    
    % Calculate and plot eigenvalues of A
    if calc_eigenvalues
%         e = eig(full(A));     % Only works for very low res grids
        e = eigs(A);    % Returns the 6 highest magnitude eigenvalues of A
        real_e = real(e);
        imag_e = imag(e);
        disp(['Max real eigenvalue: ', num2str(max(real_e))])
        disp(['Min real eigenvalue: ', num2str(min(real_e))])
        disp(['Max imag eigenvalue: ', num2str(max(abs(imag_e)))])
        disp(['Min imag eigenvalue: ', num2str(min(abs(imag_e)))])
        disp(['Max magnitude eigenvalue: ', num2str(max(abs(e)))])
        disp(['Min magnitude eigenvalue: ', num2str(min(abs(e)))])
        disp('Eigenvalues Done');
        plot(real_e, imag_e, 'ko');
        xlim([-max(abs(e)) max(abs(e))])
        ylim([-max(abs(e)) max(abs(e))])
    end
    
    % Set initial values
    u = zeros(2*m, 1);      % Initial pressure deviation
    t = 0;                  % Initial time
    
    % =====================================================================
    % INFOSTRING
    disp(['Frequency: ', num2str(w/(2*pi)), ' Hz']);
    disp(['Number of gridpoints: ', num2str(size(X_vec))])
    disp(['Simulation time: ', num2str(T), 's'])
    disp(['Number of steps: ', num2str(m_t)])
    
    % Generate id for this test
    key = join(string(randi(9,4,1)));
    key = strrep(key,' ','');
    infostring = string(append(key, '__', num2str(f), 'Hz_', num2str(m),...
        'points_', num2str(m_t), 'steps_'));
    disp(append('Test: ', num2str(f), 'Hz_', key));
    
    % Create location for this test
    location = append('../Testdata/', num2str(f), 'Hz_', num2str(key));
    
    % Saving all general data regarding this test
    if save_time_steps
        mkdir(location);
        sim_name = append(location, '/INFO.mat');
        save(sim_name, 'key', 'f', 'X_vec', 'Y_vec', 'Z_vec', 'h_t', ...
            'm_t', 'm_x', 'm_y', 'm_z', 'm', 'L_x', 'L_y', 'L_z', ...
            'infostring')
    end
    
    % =====================================================================
    % Step through time and plot
    
    % Initialize 1D plot
    if plot_time_steps_1d
        % Index to plot
        zlayer = 0.5;
        ylayer = 0.25;       % Which y-layer to plot
        lower = (round(m*zlayer, 0)-m_x*round(ylayer*m_y, 0))+...
            round(0.5*m_x*m_y);
        upper = (round(m*zlayer, 0)-m_x*round(ylayer*m_y, 0))+...
            round(0.5*m_x*m_y)+m_x-1;
        
        % Reference propagation
        propagation = amplitude(x_vec, Pvol, R);
        
        % Simple line plot with solution and references
        figure('Name', 'Pressure time plot 1D');
        u_plot = u(lower:upper);
        plot(x_vec, u_plot, 'r');
        hold on;
        plot(x_vec, propagation, 'k');
        plot(x_vec, -propagation, 'k');
        xline(1);
        hold off;
        axis([x_l x_r z]);
        title('Time: 0.0000 s');
        ylabel('Sound Pressure (Pa)');
        pause(1);
    end 
    
    % Initialize 2D plot
    if plot_time_steps_2d
        % Create 2D plot
        figure('Name', 'Pressure time plot');
        [X_vec_plot, Y_vec_plot] = meshgrid(x_vec, y_vec);
        srf = surf(X_vec_plot, Y_vec_plot, reshape(u(round(0.5*m_z,0)*...
            m_x*m_y:(round(0.5*m_z,0)+1)*m_x*m_y-1), m_y, m_x));
        axis([x_l x_r y_l y_r z]);
        pbaspect([L_x L_y min([L_x, L_y])]);
        title('Time: 0.0000 s');
        zlabel('Sound Pressure (Pa)');
        
        % Add colorbar
        cb = colorbar;
        caxis(z);
        cb.Label.String = 'Sound Pressure (Pa)';
        pause(1);
    end 
    
    % Step through time with RK4
    for time_step = 1:m_t
        [u,t] = steprk4(u, t, h_t);
        
        % Alert every 100th timestep
        if mod(time_step, 100) == 0
            disp(time_step)
        end
        
        % Save time steps
        if save_time_steps
            p = reshape(u(1:m), m_x, m_y, m_z);
            stepname = append(location, '/', num2str(key), '_', ...
                num2str(time_step), '.mat');
            save(stepname, 'p');
        end
        
        % Update 1D plot every s time steps
        if plot_time_steps_1d && mod(time_step, s) == 0
            u_plot = u(lower:upper);
            plot(x_vec, u_plot, 'r');
            hold on
            plot(x_vec, propagation, 'k');
            plot(x_vec, -propagation, 'k');
            xline(1);
            hold off
            axis([x_l x_r z]);
            ylabel('Sound Pressure (Pa)');
            title(['Time: ', num2str(time_step*h_t, '%05.4f'), ' s']);
            drawnow;
        end
        
        % Update 2D plot every s time steps
        if plot_time_steps_2d && mod(time_step, s) == 0
            % Plot middle layer
            srf.ZData = transpose(reshape(u((round(0.5*m_z,0))*m_x*m_y+1:...
                (round(0.5*m_z,0)+1)*m_x*m_y), m_x, m_y));
            srf.CData = transpose(reshape(u((round(0.5*m_z,0))*m_x*m_y+1:...
                (round(0.5*m_z,0)+1)*m_x*m_y), m_x, m_y));
            % Plot bottom layer
%             srf.ZData = transpose(reshape(u((round(1*m_z,0)-1)*m_x*m_y+1:...
%                 (round(1*m_z,0))*m_x*m_y), m_x, m_y));
%             srf.CData = transpose(reshape(u((round(1*m_z,0)-1)*m_x*m_y+1:...
%                 (round(1*m_z,0))*m_x*m_y), m_x, m_y));
            title(['Time: ', num2str(time_step*h_t, '%05.4f'), ' s']);
            drawnow;
        end
    end
    
    disp(u(1:5));
    u = reshape(u(1:m), m_x, m_y, m_z);
    disp(u(1:5,1,1));
    u = reshape(u, numel(u), 1);
    disp(u(1:5))
    
    
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





