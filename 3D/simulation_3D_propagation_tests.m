% 3D simulation of the acoustic wave equation
% 
% Simulation of the acoustic wave equation in
% three dimensions. Plots sound pressure over 
% time. 

% TO SAVE YOUR RESULTS
%   - Create a folder one step outside acoustics-kand named Testdata
%   - Your tests will be saved there (if true) in a subfolder called
%     frequencyHz_key, where key is a random 4 digit nr

% ====================================================

function simulation_3D_propagation_tests()
    
    plot_time_steps_1d = true;     % If true, plot time-steps
    plot_time_steps_2d = false;
    save_time_steps = false;    % If true, save time-steps
    
    % ====================================================
    % Model parameters
    T = 0.1;           % Final time (seconds)
    s = 1;           % plot every s time-steps
    
    % Define boundaries (m)
    x_l = 0;           % Left boundary of x
    x_r = 3;            % Right boundary of x
    L_x = x_r-x_l;      % Length of x interval
    y_l = 0;           % Left boundary of y
    y_r = 2;            % Right boundary of y
    L_y = y_r-y_l;      % Length of y interval
    z_l = 0;           % Left boundary of z
    z_r = 2;
    %z_r = 2;            % Right boundary of z
    L_z = z_r-z_l;      % Length of z interval

    % Number of grid points
    % m_x = 185;
    % m_y = 113;
    % m_z = 71;
    m_x = 51;
    m_y = 51;
    m_z = 51;
    m = m_x*m_y*m_z;

    % ====================================================
    % PDE parameters

    c = 343;              % Wave speed (m/s)
    beta_2 = c;
    beta_3 = 1;            % Absorption
    beta_3_ps = 1;

    % ====================================================
    % Initial condition parameters
    
    %lambda = max([L_x L_y L_z]); % Longest wave resonant with the room
    %k = 2;                      % Which overtone
    %f = k*c/lambda;              % Frequency
    f = 300;
    w = 2*pi*f;               % Angular frequency (room resonance)
    %w = 2*w/1.73;                 % Angular frequency (no resonance)
    %w = w*1.125;
    
    % ====================================================
    % SBP-SAT approximation

    % Spatial discretization
    h_x = L_x / (m_x-1);
    x_vec = linspace(x_l, x_r, m_x);
    h_y = L_y / (m_y-1);
    y_vec = linspace(y_l, y_r, m_y);
    h_z = L_z / (m_z-1);
    z_vec = linspace(z_l, z_r, m_z);
    [X_vec, Y_vec, Z_vec] = meshgrid(x_vec, y_vec, z_vec);
    
    % Time discretization
    h_t = 0.15*max([h_x, h_y, h_z])/c;
    m_t = round(T/h_t,0);
    h_t = T/m_t;
    
    % Amplitude of point sources
    %amp = 4*pi*34;                   % Amplitude
    %amp = 4*pi*15;
    %amp = 1/(50*h_x*h_y*h_z);     % Amplitude right at f=200, m=51, h_t=0.25
    %amp = 2/(3*h_x*h_z);            % Löjligt bra at f=200, m=51, h_t=0.25
    %amp = 1/(100000*h_x*h_z*h_t);  % bra at f=200, m=31,51,71, h_t=0.1
    amp = 3/(4*pi*h_y*h_z);
    amp_ps = amp;

    disp('Discretization Done')

    % Get D2 operator - x
    [~, HI_x, ~, D2_x, e_lx, e_rx, d1_lx, d1_rx] = sbp_cent_4th(m_x, h_x);
    % SBP-SAT
    D_x = c^2*(D2_x + HI_x*e_lx'*d1_lx - HI_x*e_rx'*d1_rx);
    E_x = -c^2/beta_2*HI_x*(beta_3_ps*e_lx'*e_lx + beta_3*e_rx'*e_rx);

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
    
%     element = sparse(m_y, m_y);
%     element(round(m_y/2), round(m_y/2)) = 1;
%     E_test = sparse(kron(element, E_x));
%     element = sparse(m_z, m_z);
%     element(round(m_z/2), round(m_z/2)) = 1;
%     E_test = sparse(kron(element, E_test));
%     
%     E = E - E_test;
    
    disp('E-Operator Done')
    
    % Construct matrix A: u_t = Au with u = [phi, phi_t]^T
    % [0, I;
    %  D, E]
    A = sparse(2*m,2*m);
    A(1:m, m+1:end) = speye(m);
    A(m+1:end, 1:m) = D;
    A(m+1:end, m+1:end) = E;
    disp('A-Matrix Done')
%     e = eig(full(E));
%     real_e = min(real(e));
%     imag_e = min(abs(imag(e)));
%     disp('Eigenvalues found');
%     plot(real_e, imag_e, 'ko');
    
    % Set initial values
    [X_vec_plot, Y_vec_plot] = meshgrid(x_vec, y_vec);
    u = zeros(2*m, 1);
    %u((round(m_x*m_y*m_z/2, 0)-m_x*round(0.5*m_y, 0))+round(0.5*m_x*m_y)) = amp;
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
    disp(append('Test: ', num2str(f), 'Hz_', key));
    
    % Create folder for this test
    location = append('../Testdata/', num2str(f), 'Hz_', num2str(key));
    
%     x_mesh = reshape(X_vec(round(0.5*m_z,0)*m_x*m_y+1:(round(0.5*m_z,0)+1)*m_x*m_y), m_y, m_x);
%     y_mesh = reshape(Y_vec(round(0.5*m_z,0)*m_x*m_y+1:(round(0.5*m_z,0)+1)*m_x*m_y), m_y, m_x);
%     %x_mesh = x_mesh(round(0.5*m_y),:);
%     %y_mesh = y_mesh(round(0.5*m_y),:);
%     x_mesh = x_mesh(:,round(0.5*m_x));
%     y_mesh = y_mesh(:,round(0.5*m_x));
%     disp(x_mesh)
%     disp(y_mesh)
    
    % Saving all general data regarding this test
    if save_time_steps
        mkdir(location);
        sim_name = append(location, '/INFO.mat');
        save(sim_name, 'key', 'f', 'X_vec', 'Y_vec', 'Z_vec', 'h_t', 'm_t', 'm_x', 'm_y', 'm_z', 'm', 'L_x', 'L_y', 'L_z', 'infostring')
    end
    
    % ====================================================
    % Plot and time step
    lower = (round(m_x*m_y*m_z/2, 0)-m_x*round(0.5*m_y, 0))+round(0.5*m_x*m_y);
    upper = (round(m_x*m_y*m_z/2, 0)-m_x*round(0.5*m_y, 0))+round(0.5*m_x*m_y)+m_x-1;
    propagation = amplitude(x_vec);
    
    % Initialize plots
    if plot_time_steps_1d
        figure('Name', 'Pressure time plot');
%         u_plot = reshape(u(round(0.5*m_z,0)*m_x*m_y:(round(0.5*m_z,0)+1)*m_x*m_y-1), m_y, m_x);
%         u_plot = u_plot(round(0.5*m_y),:);
%         %u_plot = u_plot(:,round(0.5*m_x));
%         plot(x_vec, u_plot);
%         %plot(y_vec, u_plot);
        u_plot = u(lower:upper);
        plot(x_vec, u_plot, 'r');
        hold on;
        plot(x_vec, propagation, 'k');
        plot(x_vec, -propagation, 'k');
        xline(1);
        hold off;
        z = [-10 10];
        axis([x_l x_r z]);
        %pbaspect([L_x L_y min([L_x, L_y])]);

%         % Create x and y over the slicing plane
%         xq=linspace(0,2,100);
%         yq=linspace(1,1,100);
%         disp(xq);
%         disp(yq);
% 
%         % Interpolate over the surface
%         zq=interp2(X_vec_plot,Y_vec_plot,u_plot,xq,yq); 
%         dq=sqrt((xq-0).^2 + (yq-15).^2);
% 
%         plot(dq,zq)
% 
%         axis([min(dq),max(dq),-1,1]) % to mantain a good perspective
        
        title('Time: 0 s');
        ylabel('Sound Pressure');
        
        % Add colorbar
        % cb = colorbar;
        % caxis([-1.5, 1.5]);
        % cb.Label.String = 'Sound Pressure';
        pause(1);
    end 
    % Initialize plot
    if plot_time_steps_2d
        figure('Name', 'Pressure time plot');
        srf = surf(X_vec_plot, Y_vec_plot, reshape(u(round(0.5*m_z,0)*m_x*m_y:(round(0.5*m_z,0)+1)*m_x*m_y-1), m_y, m_x));
        z = [-1 1];
        axis([x_l x_r y_l y_r z]);
        pbaspect([L_x L_y min([L_x, L_y])]);
        title('Time: 0 s');
        zlabel('Sound Pressure');
        
        % Add colorbar
        cb = colorbar;
        caxis([-1, 1]);
        cb.Label.String = 'Sound Pressure';
        pause(1);
    end 
    
    % Step through time with RK4
    for time_step = 1:m_t
        [u,t] = steprk4(u, t, h_t);
        %u = F2(t, u);
        
        if save_time_steps
            p = reshape(u(1:m), m_x, m_y, m_z);
            stepname = append(location, '/', num2str(key), '_', num2str(time_step), '.mat');
            save(stepname, 'p');
        end
        
        % Alert every 100th timestep
        if mod(time_step, 100) == 0
            disp(time_step)
        end
        
        % Plot every *insert number* time steps
        if plot_time_steps_1d && mod(time_step, s) == 0
%             u_plot = reshape(u(round(0.5*m_z,0)*m_x*m_y:(round(0.5*m_z,0)+1)*m_x*m_y-1), m_y, m_x);
%             u_plot = u_plot(round(0.5*m_y),:);
%             %u_plot = u_plot(:,round(0.5*m_x));
%             plot(x_vec, u_plot);
%             %plot(y_vec, u_plot);
            u_plot = u(lower:upper);
            plot(x_vec, u_plot, 'r');
            hold on
            plot(x_vec, propagation, 'k');
            plot(x_vec, -propagation, 'k');
            xline(1);
            hold off
            z = [-10 10];
            axis([x_l x_r z]);
            ylabel('Sound Pressure');
            title(['Time: ', num2str(time_step*h_t, '%05.4f'), ' s']);
            %pause(0.1)
            drawnow;
        end
        % Plot every *insert number* time steps
        if plot_time_steps_2d && mod(time_step, s) == 0
            % Plot middle layer
            srf.ZData = transpose(reshape(u((round(0.5*m_z,0))*m_x*m_y+1:(round(0.5*m_z,0)+1)*m_x*m_y), m_x, m_y));
            srf.CData = transpose(reshape(u((round(0.5*m_z,0))*m_x*m_y+1:(round(0.5*m_z,0)+1)*m_x*m_y), m_x, m_y));
            % Plot bottom layer
%             srf.ZData = transpose(reshape(u((round(1*m_z,0)-1)*m_x*m_y+1:(round(1*m_z,0))*m_x*m_y), m_x, m_y));
%             srf.CData = transpose(reshape(u((round(1*m_z,0)-1)*m_x*m_y+1:(round(1*m_z,0))*m_x*m_y), m_x, m_y));
            title(['Time: ', num2str(time_step*h_t, '%05.4f'), ' s']);
            drawnow;
        end
    end
    
    
    % ====================================================
    % Define functions used in code 
    
    function a = amplitude(r)
        a = 1./r;
    end

    % Define rhs of the semi-discrete approximation
    function u_t = rhs(t, u) 
        u_t = A*F2(t, u);
    end

    function v = F2(t, v)
        if t < 0.01
            amp_ps = amp*t/0.01;
        end
        v((round(m_x*m_y*m_z/2, 0)-m_x*round(0.5*m_y, 0))+round(0.5*m_x*m_y)) = amp_ps*sin(w*t);
    end

    % Time step with rk4
    function [v, t] = steprk4(v, t, dt)
%         v = F2(t, v);
        k1 = dt*rhs(t, v);
%         v = F2(t+dt/4, v);
        k2 = dt*rhs(t+0.5*dt, v+0.5*k1);
%         v = F2(t+dt/2, v);
        k3 = dt*rhs(t+0.5*dt, v+0.5*k2);
%         v = F2(t+3*dt/4, v);
        k4 = dt*rhs(t+dt, v+k3);
%         v = F2(t+dt, v);

%         k1 = dt*rhs(F2(t, v));
%         k2 = dt*rhs(F2(t+0.5*dt, v+0.5*k1));
%         k3 = dt*rhs(F2(t+0.5*dt, v+0.5*k2));
%         k4 = dt*rhs(F2(t

        v = v + 1/6*(k1 + 2*k2 + 2*k3 + k4);
        t = t + dt;
        v = F2(t, v);
    end
end