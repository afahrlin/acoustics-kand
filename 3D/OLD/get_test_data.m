

function get_test_data()
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
    
    % Time discretization
    T = 2;              % End time
    h_t = 0.01;         % Step size
    m_t = T/h_t;        % Number of steps
    
    % Define solution matrix (contains all time steps)
    u = zeros(m_z, m_x, m_y, m_t);
    t = 0;
    
    % Time step, save all values
    for i = 1:m_t
        u(:,:,:,i) = f(X_vec, Y_vec, Z_vec, t);
        t = t + h_t;
        disp(i)
    end
    
    save('test.mat', 'X_vec', 'Y_vec', 'Z_vec', 'u', 'h_t', 'm_t', 'L_x', 'L_y', 'L_z', "-v7.3")
    
    
    
    function u = g(t)
        u = sin(10*t);
    end

    function u = f(X, Y, Z, t)
        r = sqrt((X-x_0).^2+(Y-y_0).^2+(Z-z_0).^2);
        u = g(t-r);
    end
end









