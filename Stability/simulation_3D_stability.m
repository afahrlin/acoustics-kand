% Copy of final code to get eigenvalues and investigate stability
% =========================================================================

function e = simulation_3D_stability(dim, R)
    calc_eigenvalues = true;
    plot_time_steps_1d = false;     % If true, plot time steps 1d
    plot_time_steps_2d = false;     % If true, plot time steps 2d
    save_time_steps = false;        % If true, save time steps
    
    % =====================================================================
    % Model parameters
    
    T = 0;            % Final time (seconds)
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
    m_x = dim(1);
    m_y = dim(2);
    m_z = dim(3);
    m = m_x*m_y*m_z;
    
    % =====================================================================
    % PDE parameters
    
    c = 343;            % Wave speed (m/s)
    beta_2 = c;
    
    % Determine beta_3 based on reflection coefficient
%     R = 0;              % Reflection rate, must be between 0 and 1
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

    % Time discretization
%     h_t = 0.1*max([h_x, h_y, h_z])/c;
    h_t = 0.001*max([h_x, h_y, h_z]);
    m_t = round(T/h_t,0);
    h_t = T/m_t;

    % =====================================================================
    % Point source parameters
    
    % Point source frequency
%     lambda = max([L_x L_y L_z]);    % Longest wave resonant with the room
%     k = 3;                          % Which overtone
%     f = k*c/lambda;                 % Frequency
    f = 200;        % Frequency
    w = 2*pi*f;     % Angular frequency
    
    % Point source amplitude
    vol = 80;                           % Speaker volume
    Pvol = 20*10^(-6) * 10^(vol/20);    % Pressure 1m from speaker
    K = 17.0852;
    amp = Pvol/(K*h_y*h_z);        % Amplitude of point source
    amp_ps = amp;                       % Constant amp of point source
    
    % Plot amplitude
    z = [-Pvol Pvol];

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
    
    E_xy = sparse(kron(speye(m_y), E_x) + kron(E_y, speye(m_x)));
    E = sparse(kron(speye(m_z), E_xy) + kron(E_z, speye(m_x*m_y)));
    
    % Construct matrix A: v_t = Av with v = [p, p_t]^T
    % [0, I;
    %  D, E]
    A = [sparse(m,m), speye(m); D, E];
    
    % Plot eigenvalues
    if calc_eigenvalues
%         e = eig(full(A));
        e = eigs(A);
%         real_e = real(e);
%         imag_e = imag(e);
%         disp(['Max real eigenvalue: ', num2str(max(real_e))])
%         disp(['Min real eigenvalue: ', num2str(min(real_e))])
%         disp(['Max imag eigenvalue: ', num2str(max(abs(imag_e)))])
%         disp(['Min imag eigenvalue: ', num2str(min(abs(imag_e)))])
%         disp(['Max magnitude eigenvalue: ', num2str(max(abs(e)))])
%         disp(['Min magnitude eigenvalue: ', num2str(min(abs(e)))])
        e = max(abs(e));
%         disp('Eigenvalues Done');
%         figure;
%         plot(real_e, imag_e, 'ko');
%         title('Eigenvalues of Matrix A')
%         xlabel('Re(\lambda)')
%         ylabel('Im(\lambda)')
%         xline(0)
%         yline(0)
%         xlim([-max(abs(e)) max(abs(e))])
%         ylim([-max(abs(e)) max(abs(e))])
%         pbaspect([1 1 1])
%         amount = 0;
%         for i = 1:length(real_e)
%             if i == 0
%                 amount = amount + 1;
%             end
%         end
%         disp(amount)
    end
    
    
    