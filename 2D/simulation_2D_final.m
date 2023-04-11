% ====================================================
% 2D simulation of the acoustic wave equation
% 
% Simulation of the acoustic wave equation in
% two dimensions. Plots sound pressure over 
% time. 
% ====================================================

function simulation_2D_final()
    % Start timer
    tic

    % ====================================================
    % Model parameters

    plot_time_steps = true;     % If true, plot time-steps
    T = 1;                      % Final time (seconds)

    % Define boundaries (m)
    x_l = -5;           % Left boundary of x
    x_r = 5;            % Right boundary of x
    L_x = x_r-x_l;      % Length of x interval
    y_l = -5/2;         % Left boundary of y
    y_r = 5/2;          % Right boundary of y
    L_y = y_r-y_l;      % Length of y interval
    B = 1;
    a = 0.8;            % Absorption

    % Number of grid points
    m_x = 200;
    m_y = 100;
    m = m_x*m_y;

    % ====================================================
    % PDE parameters

    c = 343;              % Wave speed (m/s)

    % ====================================================
    % Initial condition parameters
    
    lambda = min([L_x L_y]);    % Shortest wave length resonant with the room
    k = 2;                      % Number of overtones (?)
    f = c/lambda;               % Frequency
    w = k*2*pi*f;               % Angular frequency (room resonance)
    %w = 1.17*w;                 % Angular frequency (no resonance)
    disp(['Frequency: ', num2str(w/(2*pi))]);
    amp = 50;                   % Amplitude

    % ====================================================
    % SBP-SAT approximation

    % Spatial discretization
    h_x = L_x / (m_x - 1);
    x_vec = linspace(x_l, x_r, m_x);
    h_y = L_y / (m_y - 1);
    y_vec = linspace(y_l, y_r, m_y);
    [X_vec, Y_vec] = meshgrid(x_vec, y_vec);

    % Time discretization
    h_t = 0.1*max([h_x, h_y])/c;
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
    
    % SBP operators
    D = sparse(kron(speye(m_y), D_x) + kron(D_y, speye(m_x)));
    E = sparse(kron(speye(m_y), E_x) + kron(E_y, speye(m_x)));

    % Construct matrix A: u_t = Au with u = [phi, phi_t]^T
    % [0, I;
    %  D, E]
    A = sparse(2*m,2*m);
    A(1:m, m+1:end) = speye(m);
    A(m+1:end, 1:m) = D;
    A(m+1:end, m+1:end) = E;

    % Set initial values (u = [phi, phi_t]^T)
    u = zeros(2*m, 1);
    t = 0;
    
    % ====================================================
    % Plot and time step
    
    % Initialize plot
    if plot_time_steps
        figure('Name', 'Pressure time plot');
        srf = surf(X_vec, Y_vec, reshape(u(1:m), m_y, m_x));
        z = [-15 15];
        axis([x_l x_r y_l y_r z]);
        pbaspect([L_x L_y min([L_x, L_y])]);
        title('Time: 0 s');
        zlabel('Sound Pressure');
        pause(1);
    end
    
    % Step through time with rk4
    for time_step = 1:m_t
        [u,t] = step(u, t, h_t);
        
        % Plot every *insert number* time steps
        if plot_time_steps && mod(time_step,20) == 0
            srf.ZData = transpose(reshape(u(1:m), m_x, m_y));
            srf.CData = transpose(reshape(u(1:m), m_x, m_y));
            title(['Time: ', num2str(time_step*h_t, '%05.4f'), ' s']);
            drawnow;
            %pause(0.001);
        end
    end

    % Stop timer
    toc
    
    % ====================================================
    % Define functions used in code

    % Define rhs of the semi-discrete approximation
    function u_t = rhs(u)
        u_t = A*u - [sparse(m,1); F(t)];
    end

    function v = F(t)
        g_lx = sparse(1, m_y);
        g_lx(round(0.25*m_y,0)) = -amp*cos(w*t);
        g_lx(round(0.75*m_y,0)) = -amp*cos(w*t);
        v = reshape((c^2*HI_x*e_lx'*g_lx), m, 1);
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


