% 3D simulation of the acoustic wave equation
% 
% Simulation of the acoustic wave equation in
% three dimensions. Plots sound pressure over 
% time. 

% ALVA FIXING
% Change time stepping method and change point sources to insertion

% TO SAVE YOUR RESULTS
%   - Create a folder one step outside acoustics-kand named Testdata
%   %% not yet - Create subfolders there, named after the frequencies you are running

% ====================================================

function simulation_3D_test()
    
    plot_time_steps = true;     % If true, plot time-steps
    save_time_steps = false;
    
    % ====================================================
    % Model parameters
    
    T = 5;           % Final time (seconds)
    
    % Define boundaries (m)
    x_l = -343;           % Left boundary of x
    x_r = 343;            % Right boundary of x
    L_x = x_r-x_l;      % Length of x interval
    y_l = -343;           % Left boundary of y
    y_r = 343;            % Right boundary of y
    L_y = y_r-y_l;      % Length of y interval
    z_l = -343;           % Left boundary of z
    z_r = 343;            % Right boundary of z
    L_z = z_r-z_l;      % Length of z interval

    % Number of grid points
    m_x = 71;
    m_y = 71;
    m_z = 71;
    m = m_x*m_y*m_z;

    % ====================================================
    % PDE parameters

    c = 343;              % Wave speed (m/s)
    beta_2 = c;
    beta_3 = 0.5;            % Absorption

    % ====================================================
    % Initial condition parameters
    
%     lambda = max([L_x L_y L_z]); % Shortest wave resonant with the room
%     k = 1;                      % Which overtone
%     f = k*c/lambda;               % Frequency
%     w = 2*pi*f;               % Angular frequency (room resonance)
%     %w = 1.17*w;                 % Angular frequency (no resonance)
%     amp = 10;                   % Amplitude
    
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
    h_t = 0.1*max([h_x, h_y, h_z])/c;
    m_t = round(T/h_t,0);
    h_t = T/m_t;

    % Get D2 operator - x
    [~, HI_x, ~, D2_x, e_lx, e_rx, d1_lx, d1_rx] = sbp_cent_6th(m_x, h_x);
    % SBP-SAT
    D_x = c^2*(D2_x + HI_x*e_lx'*d1_lx - HI_x*e_rx'*d1_rx);
    E_x = -c^2*beta_3/beta_2*HI_x*(e_lx'*e_lx + e_rx'*e_rx);

    % Get D2 operator - y
    [~, HI_y, ~, D2_y, e_ly, e_ry, d1_ly, d1_ry] = sbp_cent_6th(m_y, h_y);
    % SBP-SAT
    D_y = c^2*(D2_y + HI_y*e_ly'*d1_ly - HI_y*e_ry'*d1_ry);
    E_y = -c^2*beta_3/beta_2*HI_y*(e_ly'*e_ly + e_ry'*e_ry);
    
    % Get D2 operator - z
    [~, HI_z, ~, D2_z, e_lz, e_rz, d1_lz, d1_rz] = sbp_cent_6th(m_z, h_z);
    % SBP-SAT
    D_z = c^2*(D2_z + HI_z*e_lz'*d1_lz - HI_z*e_rz'*d1_rz);
    E_z = -c^2*beta_3/beta_2*HI_z*(e_lz'*e_lz + e_rz'*e_rz);
    
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
    u = [u_0(X_vec, Y_vec, Z_vec); zeros(m, 1)];
    t = 0;
    
    % ====================================================
    % INFOSTRING
    
    % ====================================================
    % Plot and time step
    
    % Initialize plot
    if plot_time_steps 
        figure('Name', 'Pressure time plot');
        srf = surf(X_vec_plot, Y_vec_plot, reshape(u(round(0.5*m_z,0)*m_x*m_y:(round(0.5*m_z,0)+1)*m_x*m_y-1), m_y, m_x));
        z = [-1 1];
        axis([x_l x_r y_l y_r z]);
        pbaspect([L_x L_y min([L_x, L_y])]);
        title('Time: 0 s');
        zlabel('Sound Pressure');
        
        % Add colorbar
        cb = colorbar;
        caxis([-0.5,0.5]);
        cb.Label.String = 'Sound Pressure';
        pause(1);
    end 
    
    % Step through time with RK4
    for time_step = 1:m_t
        [u,t] = step(u, t, h_t);
        
        if save_time_steps
            p = reshape(u(1:m), m_y, m_x, m_z);
            stepname = append('../Testdata/', infostring, num2str(time_step), '.mat');
            save(stepname, 'p');
        end
        
%         % Alert every 100th timestep
%         if mod(time_step, 100) == 0
%             disp(time_step)
%         end
        
        % Plot every *insert number* time steps
        if plot_time_steps && mod(time_step,4) == 0
            srf.ZData = transpose(reshape(u((round(0.5*m_z,0))*m_x*m_y+1:(round(0.5*m_z,0)+1)*m_x*m_y), m_x, m_y));
            srf.CData = transpose(reshape(u((round(0.5*m_z,0))*m_x*m_y+1:(round(0.5*m_z,0)+1)*m_x*m_y), m_x, m_y));             
            title(['Time: ', num2str(time_step*h_t, '%05.4f'), ' s']);
            drawnow;
            %pause(0.1);
        end
    end
    
    % Saving all general data regarding this test
    sim_name = append('../Testdata/INFO.mat');
    save(sim_name, 'key', 'X_vec', 'Y_vec', 'Z_vec', 'h_t', 'm_t', 'm_x', 'm_y', 'm_z', 'm', 'L_x', 'L_y', 'L_z', 'infostring')
    
    % ====================================================
    % Define functions used in code 

    % Define rhs of the semi-discrete approximation
    function u_t = rhs(u) 
        u_t = A*u;  %- [sparse(m,1); F2(t)]; 
    end

    function u = u_0(X, Y, Z)
        u = 20*exp(-(X/343/4).^2./0.05^2-(Y/343/4).^2./0.05^2-(Z/343/4).^2./0.05^2);
        u = permute(u, [2 1 3]);
        u = reshape(u, numel(X), 1);
    end

    function v = F(t)
        g_lx = sparse(m_y, m_z);
        g_lx(round(0.25*m_y,0), round(m_z*0.5, 0)) = -amp*cos(w*t);
        g_lx(round(0.75*m_y,0), round(m_z*0.5, 0)) = -amp*cos(w*t);
        G_lx = reshape(g_lx, 1, m_y*m_z);
        v = reshape((c^2*HI_x*e_lx'*G_lx), m, 1);
    end 

    function v = F2(t, v)
        v((round(m_x*m_y*m_z/2, 0) + round(m_x)*round(0.25*m_y, 0))+round(0.5*m_x*m_y)) = -amp*cos(w*t);
        v((round(m_x*m_y*m_z/2, 0) + round(m_x)*round(0.75*m_y, 0))+round(0.5*m_x*m_y)) = -amp*cos(w*t);
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