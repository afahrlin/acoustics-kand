% Plot 3D in 2D
% Run in command window
% simname = frequencyHz_key (string)

% TO PLOT YOUR RESULTS
%   - Data should be located one step outside git folder, in a subfolder
%   called /Testdata/simname


function plot_layer(simname)
    % Load general data
    location = append('../Testdata/', num2str(simname), '/');
    info = append(location, 'INFO.mat');
    load(info, 'key', 'X_vec', 'Y_vec', 'h_t', 'm_t', 'm_x', 'm_y', 'm_z', 'L_x', 'L_y', 'infostring');
    disp('Load done');
    
    % Set initial values
    
    X = X_vec(:,:,round(m_z/2));
    Y = Y_vec(:,:,round(m_z/2));
    
    s = 2;  % plot every s timesteps
    
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
    for time_step = 1:s:m_t
        step_file = append(location, num2str(key), '_', num2str(time_step), '.mat');

        load(step_file, 'p');
        U = reshape(p, m_y*m_x*m_z, 1);
        U = transpose(reshape(U((round(0.5*m_z,0))*m_x*m_y+1:(round(0.5*m_z,0)+1)*m_x*m_y), m_x, m_y));
        
        srf.CData = U;    
        %srf.ZData = U; 
        title(['Time: ', num2str((time_step-1)*h_t, '%05.4f'), ' s']);
        drawnow;
    end
end