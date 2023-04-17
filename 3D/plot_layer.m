% Plot 3D in 2D

% TO PLOT YOUR RESULTS
%   - Keep your data in a folder one step outside acoustics-kand, named Testdata


function plot_layer()
    load('../Testdata/INFO.mat', 'key', 'X_vec', 'Y_vec', 'Z_vec', 'U', 'h_t', 'm_t', 'm_x', 'm_y', 'm_z', 'L_x', 'L_y', 'L_z', 'infostring');
    
    % Set initial values
    disp(size(X_vec))
    % u = zeros(2*m, 1);
    
    X = X_vec(:,:,round(m_z/2));
    Y = Y_vec(:,:,round(m_z/2));
    
    % Prepare plot
    figure('Name', 'Pressure time plot');
    srf = surf(X, Y, zeros(size(X)));
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
        P = reshape(p, m_y*m_x*m_z, 1);
        P = transpose(reshape(P((round(0.5*m_z,0))*m_x*m_y+1:(round(0.5*m_z,0)+1)*m_x*m_y), m_x, m_y));
        %p = permute(p,[2,1,3]);
        %P = p(:,:,round(m_z/2));
        %P = reshape(P, numel(P), 1);

        if time_step == 1
            size(p)
            size(P)
        end

        if mod(time_step,2) == 0
            srf.ZData = P;
            srf.CData = P;             
            title(['Time: ', num2str(time_step*h_t, '%05.4f'), ' s']);
            drawnow;
        end
    end
end