% ====================================================
% 2D simulation of the acoustic wave equation
% 
% Simulation of the acoustic wave equation in
% two dimensions. Plots sound pressure over 
% time. 
% ====================================================

function simulation_2D_transposed()
    % Start timer
    tic


    % ====================================================
    % Model parameters

    plot_time_steps = true;     % If true, plot time-steps
    T = 3;                      % Final time

    % Define boundaries
    x_l = -1;%-1/2;         % Left boundary of x
    x_r = 1;%1/2;          % Right boundary of x
    L_x = x_r-x_l;      % Length of x interval
    y_l = -1/2;         % Left boundary of y
    y_r = 1/2;          % Right boundary of y
    L_y = y_r-y_l;      % Length of y interval

    % Number of grid points
    m_x = 100;
    m_y = 100;
    m = m_x*m_y;

    % ====================================================
    % PDE parameters

    c = 1;              % Wave speed

    % ====================================================
    % Initial condition parameters

    x_0 = 0;            % Center - x
    y_0 = 0;            % Center - y
    o = 0.05;           % Sigma (initial function)

    % ====================================================
    % SBP-SAT approximation

    % Spatial discretization
    h_x = L_x / (m_x - 1);
    x_vec = linspace(x_l, x_r, m_x);
    h_y = L_y / (m_y - 1);
    y_vec = linspace(y_l, y_r, m_y);
    [X_vec, Y_vec] = meshgrid(x_vec, y_vec);

    % Time discretization
    h_t = 0.1*min([h_x, h_y])/c;
    m_t = round(T/h_t,0);
    h_t = T/m_t;
    disp('Discretizations done');

    % Get D2 operator - x
    [~, HI_x, ~, D2_x, e_lx, e_rx, d1_lx, d1_rx] = sbp_cent_6th(m_x, h_x);
    % SBP-SAT
    D_x = c^2*D2_x + c^2*HI_x*transpose(e_lx)*d1_lx*eye(m_x) - c^2*HI_x*transpose(e_rx)*d1_rx*eye(m_x);

    % Get D2 operator - y
    [~, HI_y, ~, D2_y, e_ly, e_ry, d1_ly, d1_ry] = sbp_cent_6th(m_y, h_y);
    % SBP-SAT
    D_y = c^2*D2_y + c^2*HI_y*transpose(e_ly)*d1_ly*eye(m_y) - c^2*HI_y*transpose(e_ry)*d1_ry*eye(m_y);
    
    % SBP operator
    D = sparse(kron(eye(m_y), D_x) + kron(D_y, eye(m_x)));
    disp('D done');

    % Construct matrix: u_t = Au with u = [phi, phi_t]^T
    % [0, I;
    %  D, 0]
    A = zeros(2*m);
    disp('A step 1 done')
    A(1:m, m+1:end) = sparse(eye(m));
    disp('A step 2 done')
    A(m+1:end, 1:m) = D;
    disp('A step 3 done')
    A = sparse(A);
    disp('A done');

    % Set initial values (u = [phi, phi_t]^T)
    u = [phi_0(X_vec, Y_vec, m, o, x_0, y_0); zeros(1, m)];
    t = 0;
    
    %disp(size(u))
    disp(u(9995: 10005))
    
    surf(X_vec, Y_vec, reshape(u(m+1:end), m_y, m_x));
    z = [-15 15];
    axis([x_l x_r y_l y_r z]);
    pbaspect([L_x L_y min([L_x, L_y])])
    drawnow;
    pause(1);
    
    for time_step = 1:m_t
        [u,t] = step(u, t, h_t);
        if mod(time_step,10) == 0
            surf(X_vec, Y_vec, reshape(u(m+1:end), m_y, m_x));
            z = [-15 15];
            axis([x_l x_r y_l y_r z]);
            pbaspect([L_x L_y min([L_x, L_y])])
            drawnow;
            pause(0.01);
        end
    end

    % Stop timer
    toc

    % Define rhs of the semi-discrete approximation
    function u_t = rhs(u)
        u_t = A*u;
    end

    % Define initial function
    function u = phi_0(x, y, m, o, x_0, y_0)
        u = reshape(exp(-(((x-x_0).^2)./(o^2))-(((y-y_0).^2)./(o^2))), 1, m);
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























