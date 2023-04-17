% Plot 3D in 2D

% TO PLOT YOUR RESULTS
%   - Keep your data in a folder one step outside acoustics-kand, named Testdata


function plot_layer()
    load('../Testdata/INFO.mat', 'key', 'x_vec', 'y_vec', 'z_vec', 'h_t', 'm_t', 'm_x', 'm_y', 'm_z', 'm', 'L_x', 'L_y', 'L_z', 'infostring');
    
    % Set initial values
    [X_vec_plot, Y_vec_plot] = meshgrid(x_vec, y_vec);
    disp('Xvecplot')
    size(X_vec_plot)
    u = zeros(2*m, 1);
    
    % Prepare plot
    figure('Name', 'Pressure time plot');
    srf = surf(X_vec_plot, Y_vec_plot, reshape(u(round(0.5*m_z,0)*m_x*m_y:(round(0.5*m_z,0)+1)*m_x*m_y-1), m_y, m_x));
    z = [-10 10];
    axis([-L_x/2 L_x/2 -L_y/2 L_y/2 z]);
    pbaspect([L_x L_y min([L_x, L_y])]);
    title('Time: 0 s');
    zlabel('Sound Pressure');
    pause(1);
    
    % Time step
    for time_step = 1:m_t
    step_file = append('../Testdata/', infostring, num2str(time_step), '.mat');
    
    load(step_file, 'p');
    P = permute(p,[2,1,3]);
    
    if time_step == 1
        size(p)
        size(P)
    end
    
        if mod(time_step,2) == 0
            u = (P(:,:,round(m_z/2)));
            size(u)
            srf.ZData = (u);
            srf.CData = (u);             
            title(['Time: ', num2str(time_step*h_t, '%05.4f'), ' s']);
            drawnow;
        end
    end
    
end