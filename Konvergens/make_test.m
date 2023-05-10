% Making convergence test data

function [u, simname] = make_test(f, T, dim)

    save_time_steps = false;    % If true, save time-steps
    
    % ====================================================
    % Model parameters
    
    x_l = 0;           % Left boundary of x
    x_r = 6.24;            % Right boundary of x
    L_x = x_r-x_l;      % Length of x interval
    y_l = 0;           % Left boundary of y
    y_r = 3.84;            % Right boundary of y
    L_y = y_r-y_l;      % Length of y interval
    z_l = 0;           % Left boundary of z
    z_r = 2.4;            % Right boundary of z
    L_z = z_r-z_l;      % Length of z interval

    % Number of grid points
    m_x = dim(1);
    m_y = dim(2);
    m_z = dim(3);
    m = m_x*m_y*m_z;

    % ====================================================
    % PDE parameters

    c = 343;              % Wave speed (m/s)
    beta_2 = c;
    beta_3 = 0;            % Absorption

    % ====================================================
    % Initial condition parameters
    
    vol = 80;
    Pvol = 20*10^(-6) * 10^(20/vol);
    w = 2*pi*f;               % Angular frequency 
    amp = 4*pi*3*Pvol;        % Amplitude
    amp_ps = amp;
    
    % ====================================================
    % SBP-SAT approximation

    % Spatial discretization
    h_x = L_x / (m_x - 1);
    x_vec = linspace(x_l, x_r, m_x);
    h_y = L_y / (m_y - 1);
    y_vec = linspace(y_l, y_r, m_y);
    h_z = L_z / (m_z - 1);
    z_vec = linspace(z_l, z_r, m_z);
    [X_vec, Y_vec, Z_vec] = meshgrid(x_vec, y_vec, z_vec);

    % Time discretization
    h_t = 0.25*max([h_x, h_y, h_z])/c;
    m_t = round(T/h_t,0);
    h_t = T/m_t;

    disp('Discretization Done')

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

    % Construct matrix A: u_t = Au with u = [phi, phi_t]^T
    % [0, I;
    %  D, E]
    A = sparse(2*m,2*m);
    A(1:m, m+1:end) = speye(m);
    A(m+1:end, 1:m) = D;
    A(m+1:end, m+1:end) = E;
    disp('A-Matrix Done')
    
    % Set initial values
    u = zeros(2*m, 1) + 20*10^(-6);
    t = 0;
    
    % ====================================================
    % INFOSTRING
    disp(['Frequency: ', num2str(w/(2*pi)), ' Hz']);
    disp(['Number of gridpoints: ', num2str(size(X_vec))])
    disp(['Simulation time: ', num2str(T), 's'])
    disp(['Number of steps: ', num2str(m_t)])
    
    % Generate id for this test
    key = join(string(randi(9,4,1)));
    key = strrep(key,' ','');
    infostring = string(append(key, '__', num2str(f), 'Hz_', num2str(m), 'points_', num2str(m_t), 'steps_'));
    simname = append(num2str(f), 'Hz_', num2str(m), 'points');
    disp(append('Test: ', simname));
    
    % Create folder for this test
    location = append('Testdata/', simname);
    mkdir(location)

    sim_info = append(location, '/INFO.mat');
    save(sim_info, 'simname', 'key', 'f', 'X_vec', 'Y_vec', 'Z_vec', 'h_t', 'm_t', 'm_x', 'm_y', 'm_z', 'm', 'L_x', 'L_y', 'L_z', 'infostring', '-v7.3')
    disp('Parameters saved');
    % ====================================================
    
    % Step through time with RK4
    for time_step = 1:m_t
        [u,t] = steprk4(u, t, h_t);
        
        if save_time_steps
            p = reshape(u(1:m), m_x, m_y, m_z);
            stepname = append(location, '/', num2str(key), '_', num2str(time_step), '.mat');
            save(stepname, 'p');
        end
        
        % Alert every 100th timestep
        if mod(time_step, 100) == 0
            disp([num2str(time_step), '/', num2str(m_t)])
        end
    end
    
    u = reshape(u(1:m), m_x, m_y, m_z);
    
    
    % ====================================================
    % Define functions used in code 

    % Define rhs of the semi-discrete approximation
    function u_t = rhs(u) 
        u_t = A*u;
    end

    function v = F(t)
        g_lx = sparse(m_y, m_z);
        g_lx(round(0.25*m_y,0), round(m_z*0.5, 0)) = -amp*cos(w*t);
        g_lx(round(0.75*m_y,0), round(m_z*0.5, 0)) = -amp*cos(w*t);
        G_lx = reshape(g_lx, 1, m_y*m_z);
        v = reshape((c^2*HI_x*e_lx'*G_lx), m, 1);
    end 

    function v = F2(t, v)
        if t < 0.05
            amp_ps = amp*t/0.05;
        end

        % One point source in the middle
        %v(round(m_x*m_y*m_z/2, 0)+m_x*m_y) = amp*sin(w*t);
        
        % Two point sources on boundary
        v((round(m_x*m_y*m_z/2, 0)-m_x*round(0.25*m_y, 0))+round(0.5*m_x*m_y)) = amp_ps*sin(w*t);
        v((round(m_x*m_y*m_z/2, 0)-m_x*round(0.75*m_y, 0))+round(0.5*m_x*m_y)) = amp_ps*sin(w*t);
    end

    % Time step with rk4
    function [v, t] = steprk4(v, t, dt)
        v = F2(t, v);
        k1 = dt*rhs(v);
        v = F2(t+dt/4, v);
        k2 = dt*rhs(v+0.5*k1);
        v = F2(t+dt/2, v);
        k3 = dt*rhs(v+0.5*k2);
        v = F2(t+3*dt/4, v);
        k4 = dt*rhs(v+k3);

        v = v + 1/6*(k1 + 2*k2 + 2*k3 + k4);
        t = t + dt;
    end
end