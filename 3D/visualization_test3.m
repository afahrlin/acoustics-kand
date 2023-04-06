% ====================================================
% Test of visualization of some 3D data with motion
% ====================================================


function visualization_test3()
    filename = 'test.mat';
    load(filename, 'X_vec', 'Y_vec', 'Z_vec', 'u', 'h_t', 'm_t', 'L_x', 'L_y', 'L_z');
    
    % Reshape and select every nth index to visualize
    n = 3;
    X = X_vec(1:n:end,1:n:end,1:n:end);
    X = reshape(X, numel(X), 1);
    Y = Y_vec(1:n:end,1:n:end,1:n:end);
    Y = reshape(Y, numel(Y), 1);
    Z = Z_vec(1:n:end,1:n:end,1:n:end);
    Z = reshape(Z, numel(Z), 1);
    
    U = u(:,:,:,1);
    U = U(1:n:end,1:n:end,1:n:end);
    U = reshape(U, numel(U), 1);
    
    % Initialize plot ===============================
    % Plot as a scatter plot
    figure;
    sc1 = scatter3(X, Y, Z, 1, U, 'filled');
    view(-31,14)
    xlabel('X')
    ylabel('Y')
    zlabel('Z')
    title('Visualization 2! Time: 0');
    pbaspect([L_x L_y L_z]);
    
    m = 1;
    
    for i = 1:m_t/m
        U = u(:,:,:,i*m);
        U = U(1:n:end,1:n:end,1:n:end);
        U = reshape(U, numel(U), 1);
        tic
        sc1.CData = U;
        toc
        title(['Visualization 2! Time: ', num2str(h_t*(i*m-1))]);
        pause(0.001);
    end
end







