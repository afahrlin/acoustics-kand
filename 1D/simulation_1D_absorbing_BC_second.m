% ====================================================
% 2D simulation of the acoustic wave equation
% 
% Simulation of the acoustic wave equation in
% two dimensions. Plots sound pressure over 
% time. 
% ====================================================

function simulation_1D_absorbing_BC()
    % Start timer
    tic

    % ====================================================
    % Model parameters

    plot_time_steps = true;     % If true, plot time-steps
    T = 3;                      % Final time

    % Define boundaries
    x_l = -343;           % Left boundary of x
    x_r = 343;            % Right boundary of x
    L_x = x_r-x_l;      % Length of x interval


    % Number of grid points
    m = 200;

    % ====================================================
    % PDE parameters

    c = 343;              % Wave speed
    B = c;
    a = 0.5;
    
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
    h_t = 0.1*h_x/c;
    m_t = round(T/h_t,0);
    h_t = T/m_t;

    % Get D2 operator - x
    [~, HI_x, ~, D2_x, e_lx, e_rx, d1_lx, d1_rx] = sbp_cent_6th(m, h_x);
    % SBP-SAT
    D_x = c^2*D2_x + c^2*HI_x*e_lx'*d1_lx - c^2*HI_x*e_rx'*d1_rx;
    E_x = - a/B*c^2*HI_x*e_lx'*e_lx - a/B*c^2*HI_x*e_rx'*e_rx;

    % Construct matrix: u_t = Au with u = [phi, phi_t]^T
    % [0, I;
    %  D, E]
    A = sparse(2*m,2*m);
    A(m+1:end, m+1:end) = sparse(E_x);
    A(1:m, m+1:end) = speye(m);
    A(m+1:end, 1:m) = sparse(D_x);

    % Set initial values (u = [u, u_t]^T)
    u = [phi_0(X_vec); zeros(m, 1)];
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
        [u,t] = step(u, t, h_t);
        
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

    % Define rhs of the semi-discrete approximation
    function u_t = rhs(u)
        u_t = A*u;
    end

    % Define initial function
    function u = phi_0(x)
        u = exp(-((((x-x_0)/343/10).^2)./(o^2)));
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























