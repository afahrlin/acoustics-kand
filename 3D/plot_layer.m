% Plot 3D data as 2D layer at height 'height'
% simname format: frequencyHz_key (string)
% height = relative height of layer; [0 , 1]

% TO PLOT YOUR RESULTS
%   - Data should be located one step outside git folder, in a subfolder
%   called /Testdata/simname


function plot_layer(simname, height)
    % Load general data
    location = append('../Testdata/', num2str(simname), '/');
    info = append(location, 'INFO.mat');
    load(info, 'key', 'f', 'X_vec', 'Y_vec', 'h_t', 'm_t', 'm_x', 'm_y', 'm_z', 'L_x', 'L_y', 'infostring');
    disp('Load done');
    disp(infostring);
    
    % Set initial values
    
    X = X_vec(:,:,round(m_z/2));
    Y = Y_vec(:,:,round(m_z/2));
    
    s = 2;  % plot every s timesteps
    % height = 0.55;   % percentage of full height
    Height = zeros(m_x*m_y, 1) + height;
    
    % Prepare plot
    figure('Name', append('Sound pressure (layer) ', num2str(f), ' Hz'));
    srf = surf(X, Y, zeros(size(X)));
    z = [0 1];
    axis([0 L_x 0 L_y z]);
    srf.ZData = transpose(reshape(Height, m_x, m_y));
    srf.CData = transpose(reshape(zeros(m_x*m_y, 1), m_x, m_y));
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
        step_file = append(location, key, '_', num2str(time_step), '.mat');

        load(step_file, 'p');
        U = reshape(p, m_y*m_x*m_z, 1);
        U = transpose(reshape(U((round(height*m_z,0))*m_x*m_y+1:(round(height*m_z,0)+1)*m_x*m_y), m_x, m_y));
        
        srf.CData = U;
        title(['Time: ', num2str((time_step-1)*h_t, '%05.4f'), ' s']);
        drawnow;
    end
end