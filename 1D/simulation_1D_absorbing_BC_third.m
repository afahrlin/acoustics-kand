% ====================================================
% 2D simulation of the acoustic wave equation
% 
% Simulation of the acoustic wave equation in
% two dimensions. Plots sound pressure over 
% time. 
% ====================================================

function simulation_1D_absorbing_BC_third()
    % Start timer
    tic

    % ====================================================
    % Model parameters

    plot_time_steps = true;     % If true, plot time-steps
    T = 2;                      % Final time

    % Define boundaries
    x_l = -343;           % Left boundary of x
    x_r = 343;            % Right boundary of x
    L_x = x_r-x_l;      % Length of x interval


    % Number of grid points
    m = 200;

    % ====================================================
    % PDE parameters

    c = 343;              % Wave speed
    %B = 9.76*c^3;
    B = c;
    a = 1;
    
    % ====================================================
    % Initial condition parameters

    x_0 = 0;            % Center - x
    o = 0.05;           % Sigma (initial function)

    % ====================================================
    % SBP-SAT approximation

    % Spatial discretization
    h_x = L_x / (m - 1);
    X_vec = linspace(x_l, x_r, m)';

    % Time discretization
    h_t = 0.01*h_x/c;
    m_t = round(T/h_t,0);
    h_t = T/m_t;

    % Get D2 operator - x
    [~, HI_x, ~, D2_x, e_lx, e_rx, d1_lx, d1_rx] = sbp_cent_6th(m, h_x);
    % SBP-SAT
    D = c^2*D2_x + c^2*HI_x*e_lx'*d1_lx - c^2*HI_x*e_rx'*d1_rx;
    E = - a/B*c^2*HI_x*e_lx'*e_lx - a/B*c^2*HI_x*e_rx'*e_rx;

    % Set initial values (u = [u, u_t]^T)
    u = phi_0(X_vec);
    u_prev = phi_0(X_vec);
    t = 0;
    
    % ====================================================
    % Plot and time step
    
    % Initialize plot
    if plot_time_steps
        %plot(X_vec, -u(m+1:end));
        plot(X_vec, u(1:m));
        %z = [-15 15];
        z = [-2 2];
        axis([x_l x_r z]);
        drawnow;
        pause(2);
    end
    
    % Step through time with rk4
    for time_step = 1:m_t
        [u,u_prev,t] = step(u, u_prev, h_t, t);
        
        % Plot every 10 time steps
        if plot_time_steps && mod(time_step,10) == 0
            %plot(X_vec, -u(m+1:end));
            plot(X_vec, u(1:m));
            %z = [-15 15];
            z = [-2 2];
            axis([x_l x_r z]);
            drawnow;
            pause(0.01);
        end
    end

    % Stop timer
    toc
    
    % ====================================================
    % Define functions used in code

    % Define initial function
    function u = phi_0(x)
        u = 2*exp(-((((x-x_0)/343).^2)./(o^2)));
    end

    function [v_new, v, t] = step(v, v_prev, dt, t)
        vh_new = dt^2*D*v + 2*v - v_prev;
        v_new = vh_new + E*(vh_new-v_prev)/(2*dt);

        t = t + dt;
    end
end























