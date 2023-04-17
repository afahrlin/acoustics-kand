% Plot 3D in 2D

% TO PLOT YOUR RESULTS
%   - Keep your data in a folder one step outside acoustics-kand, named Testdata


function plot_layer()
    load('../Testdata/INFO.mat', 'key', 'X_vec', 'Y_vec', 'Z_vec', 'h_t', 'm_t', 'm_x', 'm_y', 'm_z', 'L_x', 'L_y', 'L_z', 'infostring');
    
    % Set initial values
    disp(size(X_vec))
    
    X = X_vec(:,:,round(m_z/2));
    Y = Y_vec(:,:,round(m_z/2));
    
    % Prepare plot
    figure('Name', 'Pressure time plot');
    srf = surf(X, Y, zeros(size(X)));
    z = [-1 1];
    axis([-L_x/2 L_x/2 -L_y/2 L_y/2 z]);
    pbaspect([L_x L_y min([L_x, L_y])]);
    title('Time: 0 s');
    zlabel('Sound Pressure');
    
    % Add colorbar
    cb = colorbar;
    caxis([-0.5,0.5]);
    cb.Label.String = 'Sound Pressure';
    pause(1);
    
    % Time step
    for time_step = 1:m_t
        step_file = append('../Testdata/', infostring, num2str(time_step), '.mat');

        load(step_file, 'p');
        P = reshape(p, m_y*m_x*m_z, 1);
        P = transpose(reshape(P((round(0.5*m_z,0))*m_x*m_y+1:(round(0.5*m_z,0)+1)*m_x*m_y), m_x, m_y));

        if mod(time_step,2) == 0
            srf.CData = P;             
            title(['Time: ', num2str(time_step*h_t, '%05.4f'), ' s']);
            drawnow;
        end
    end
end