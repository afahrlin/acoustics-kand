% 3D simulation of the acoustic wave equation
% 
% Simulation of the acoustic wave equation in
% three dimensions.

% TO SAVE YOUR RESULTS
%   - Create a folder one step outside acoustics-kand named Testdata
%   - Your tests will be saved there (if true) in a subfolder called
%     frequencyHz_key, where key is a random 4 digit nr

% ====================================================

function simulation_3D_second()
    
    plot_time_steps = false;    % If true, plot time-steps
    save_time_steps = true;     % If true, save time-steps
    
    % ====================================================
    % Model parameters
    
    T = 0.05;             % Final time (seconds)
    s = 2;                % plot every s time-steps
    
    % Define boundaries (m)
    x_l = -7/2;           % Left boundary of x
    x_r = 7/2;            % Right boundary of x
    L_x = x_r-x_l;        % Length of x interval
    y_l = -5/2;           % Left boundary of y
    y_r = 5/2;            % Right boundary of y
    L_y = y_r-y_l;      % Length of y interval
    z_l = -3;           % Left boundary of z
    z_r = 3;            % Right boundary of z
    L_z = z_r-z_l;      % Length of z interval
    
    B = 1;
    a = 0.8;            % Absorption

    % Number of grid points
    m_x = 61;
    m_y = 51;
    m_z = 41;
    m = m_x*m_y*m_z;

    % ====================================================
    % PDE parameters

    c = 343;              % Wave speed (m/s)

    % ====================================================
    % Initial condition parameters
    
    lambda = min([L_x L_y L_z]); % Shortest wave resonant with the room
    k = 2;                      % Which overtone
    f = c/lambda;               % Frequency
    w = k*2*pi*f;               % Angular frequency (room resonance)
    %w = 1.17*w;                 % Angular frequency (no resonance)
    amp = 50;                   % Amplitude
    
    % ====================================================
    % SBP-SAT approximation

    % Spatial discretization
    h_x = L_x / (m_x - 1);
    x_vec = linspace(x_l, x_r, m_x);
    h_y = L_y / (m_y - 1);
    y_vec = linspace(y_l, y_r, m_y);
    h_z = L_z / (m_z - 1);
    z_vec = linspace(y_l, y_r, m_z);
    [X_vec, Y_vec, Z_vec] = meshgrid(x_vec, y_vec, z_vec);

    % Time discretization
    h_t = 0.25*max([h_x, h_y, h_z])/c;
    m_t = round(T/h_t,0);
    h_t = T/m_t;

    % Get D2 operator - x
    [~, HI_x, ~, D2_x, e_lx, e_rx, d1_lx, d1_rx] = sbp_cent_6th(m_x, h_x);
    % SBP-SAT
    D_x = c^2*D2_x + c^2/B*HI_x*e_lx'*d1_lx - c^2/B*HI_x*e_rx'*d1_rx;
    E_x = - a/B*HI_x*e_lx'*e_lx - a/B*HI_x*e_rx'*e_rx;

    % Get D2 operator - y
    [~, HI_y, ~, D2_y, e_ly, e_ry, d1_ly, d1_ry] = sbp_cent_6th(m_y, h_y);
    % SBP-SAT
    D_y = c^2*D2_y + c^2/B*HI_y*e_ly'*d1_ly - c^2/B*HI_y*e_ry'*d1_ry;
    E_y = - a/B*HI_y*e_ly'*e_ly - a/B*HI_y*e_ry'*e_ry;
    
    % Get D2 operator - z
    [~, HI_z, ~, D2_z, e_lz, e_rz, d1_lz, d1_rz] = sbp_cent_6th(m_z, h_z);
    % SBP-SAT
    D_z = c^2*D2_z + c^2/B*HI_z*e_lz'*d1_lz - c^2/B*HI_z*e_rz'*d1_rz;
    E_z = - a/B*HI_z*e_lz'*e_lz - a/B*HI_z*e_rz'*e_rz;
    
    % SBP operators
    D_xy = sparse(kron(speye(m_y), D_x) + kron(D_y, speye(m_x)));
    D = sparse(kron(speye(m_z), D_xy) + kron(D_z, speye(m_x*m_y)));
    
    E_xy = sparse(kron(speye(m_y), E_x) + kron(E_y, speye(m_x)));
    E = sparse(kron(speye(m_z), E_xy) + kron(E_z, speye(m_x*m_y)));

    % Construct matrix A: u_t = Au with u = [phi, phi_t]^T
    % [0, I;
    %  D, E]
    A = sparse(2*m,2*m);
    A(1:m, m+1:end) = speye(m);
    A(m+1:end, 1:m) = D;
    A(m+1:end, m+1:end) = E;
    
    % Set initial values
    [X_vec_plot, Y_vec_plot] = meshgrid(x_vec, y_vec);
    u = zeros(2*m, 1);
    t = 0;
    
    % ====================================================
    % INFOSTRING
    disp(['Frequency: ', num2str(f), ' Hz']);
    disp(['Number of gridpoints: ', num2str(size(X_vec))])
    disp(['Simulation time: ', num2str(T), 's'])
    disp(['Number of steps: ', num2str(m_t)])
    
    key = join(string(randi(9,4,1)));
    key = strrep(key,' ','');
    infostring = string(append(key, '__', num2str(f), 'Hz_', num2str(m), 'points_', num2str(m_t), 'steps_'));
    
    location = append('../Testdata/', num2str(f), 'Hz_', num2str(key));
    mkdir(location);   % create the directory
    
    % ====================================================
    % Plot and time step
    
    % Initialize plot
    if plot_time_steps 
        figure('Name', 'Pressure time plot');
        srf = surf(X_vec_plot, Y_vec_plot, reshape(u(round(0.5*m_z,0)*m_x*m_y:(round(0.5*m_z,0)+1)*m_x*m_y-1), m_y, m_x));
        z = [-10 10];
        axis([x_l x_r y_l y_r z]);
        pbaspect([L_x L_y min([L_x, L_y])]);
        title('Time: 0 s');
        zlabel('Sound Pressure');
        pause(1);
    end 
    
    % Step through time with RK4
    for time_step = 1:m_t
        [u,t] = step(u, t, h_t);
        
        if save_time_steps
            p = reshape(u(1:m), m_x, m_y, m_z);
            stepname = append(location, '/', num2str(key), '_', num2str(time_step), '.mat');
            save(stepname, 'p');
        end
        
        % Plot every s time steps
        if plot_time_steps && mod(time_step, s) == 0
            srf.ZData = transpose(reshape(u((round(0.5*m_z,0))*m_x*m_y+1:(round(0.5*m_z,0)+1)*m_x*m_y), m_x, m_y));
            srf.CData = transpose(reshape(u((round(0.5*m_z,0))*m_x*m_y+1:(round(0.5*m_z,0)+1)*m_x*m_y), m_x, m_y));             
            title(['Time: ', num2str(time_step*h_t, '%05.4f'), ' s']);
            drawnow;
        end
    end
    
    % Saving all general data regarding this test
    sim_name = append(location, '/INFO.mat');
    save(sim_name, 'key', 'X_vec', 'Y_vec', 'Z_vec', 'h_t', 'm_t', 'm_x', 'm_y', 'm_z', 'm', 'L_x', 'L_y', 'L_z', 'infostring')
    
    % ====================================================
    % Define functions used in code 

    % Define rhs of the semi-discrete approximation
    function u_t = rhs(u) 
        u_t = A*u  - [sparse(m,1); F(t)]; 
    end

    function v = F(t)
        g_lx = sparse(m_y, m_z);
        g_lx(round(0.25*m_y,0), round(m_z*0.5, 0)) = -amp*cos(w*t);
        g_lx(round(0.75*m_y,0), round(m_z*0.5, 0)) = -amp*cos(w*t);
        G_lx = reshape(g_lx, 1, m_y*m_z);
        v = reshape((c^2*HI_x*e_lx'*G_lx), m, 1);
    end 

    % Time step with rk4
    function [v, t] = step(v, t, dt)
        k1 = dt*rhs(v);
        k2 = dt*rhs(v+0.5*k1);
        k3 = dt*rhs(v+0.5*k2);
        k4 = dt*rhs(v+k3);

        v = v + 1/6*(k1 + 2*k2 + 2*k3 + k4);
        t = t + dt;
    end
end