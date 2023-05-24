% test convergence based on UPPMAX runs
% L2 norm

% Load everything
function conv_study_olle()
%     uref = struct2array(load('Konvergens/REF/u.mat', 'u'));
    uref = struct2array(load('../Testdata/conv_test/REF/u.mat', 'u'));
%     uref = permute(uref, [2 1 3]);
    % normalising uref (illegal?)
%     uref = uref./max(max(max(uref)));

    u1 = struct2array(load(append('../Testdata/conv_test/200Hz_36729points__1/u.mat'), 'u'));
    u2 = struct2array(load(append('../Testdata/conv_test/200Hz_279825points__2/u.mat'), 'u'));
    u3 = struct2array(load(append('../Testdata/conv_test/200Hz_928969points__3/u.mat'), 'u'));
    u4 = struct2array(load(append('../Testdata/conv_test/200Hz_2183841points__4/u.mat'), 'u'));
    u6 = struct2array(load(append('../Testdata/conv_test/200Hz_7309489points__6/u.mat'), 'u'));
%     
%     u1 = struct2array(load(append('Konvergens/t1/u.mat'), 'u'));
%     u2 = struct2array(load(append('Konvergens/t2/u.mat'), 'u'));
%     u3 = struct2array(load(append('Konvergens/t3/u.mat'), 'u'));
%     u4 = struct2array(load(append('Konvergens/t4/u.mat'), 'u'));
%     u6 = struct2array(load(append('Konvergens/t6/u.mat'), 'u'));
   
    n = [12, 6, 4, 3, 2];
    h = [0.12, 0.06, 0.04, 0.03, 0.02];
    l2_norm = zeros(1, 5);
    H_norm = zeros(1, 5);

    for i = 1:5
        if i == 1
            u=u1;
        elseif i == 2
            u=u2;
        elseif i == 3
            u=u3;
        elseif i == 4
            u=u4;
        elseif i == 5
            u=u6;
        end
        
%         u = permute(u, [2 1 3]);
        
        % normalising u, (illegal?)
%         u = u./max(max(max(u)));

        % values relevant to this size
        ref = uref(1:n(i):end, 1:n(i):end, 1:n(i):end);
%         disp(size(u))
%         disp(size(ref))

        % difference
        diff = ref - u;
        diff = reshape(diff, numel(diff), 1);
        
        m_x = size(u, 1);
        m_y = size(u, 2);
        m_z = size(u, 3);
        h_x = 6.24 / (m_x-1);
        h_y = 3.84 / (m_y-1);
        h_z = 2.4 / (m_z-1);
        [H_x, ~, ~, ~, ~, ~, ~, ~] = sbp_cent_4th(m_x, h_x);
        [H_y, ~, ~, ~, ~, ~, ~, ~] = sbp_cent_4th(m_y, h_y);
        [H_z, ~, ~, ~, ~, ~, ~, ~] = sbp_cent_4th(m_z, h_z);
%         H_x = h_x*speye(m_x);
%         H_y = h_y*speye(m_y);
%         H_z = h_z*speye(m_z);
        
%         H_xy = h(i).*sparse(kron(speye(m_y), H_x) + kron(H_y, speye(m_x)));
%         H = h(i).*sparse(kron(speye(m_z), H_xy) + kron(H_z, speye(m_x*m_y)));
%         H_xy = sparse(h_y*kron(speye(m_y), H_x) + h_x*kron(H_y, speye(m_x)));
%         H = sparse(h_z*kron(speye(m_z), H_xy) + h_x*h_y*kron(H_z, speye(m_x*m_y)));
%         H_xy = sparse(kron(speye(m_y), H_x*h_y*h_z) + kron(H_y*h_x*h_z, speye(m_x)));
%         H = sparse(kron(speye(m_z), H_xy) + kron(H_z*h_x*h_y, speye(m_x*m_y)));
%         H_xy = sparse(kron(speye(m_y), H_x/h_x) + kron(H_y/h_y, speye(m_x)))/2;
%         H = sparse(kron(speye(m_z), H_xy) + kron(H_z/h_z, speye(m_x*m_y)))/2;
        
%         H_xy = sparse(kron(h_y*speye(m_y), H_x) + kron(H_y, speye(m_x)*h_x));
%         H = sparse(kron(h_z*speye(m_z), H_xy) + kron(H_z, h_x*h_y*speye(m_x*m_y)));
%         H = (sparse(kron(kron(h_z*speye(m_z), h_y*speye(m_y)), H_x)) + ...
%             sparse(kron(kron(h_z*speye(m_z), H_y), h_x*speye(m_x))) + ...
%             sparse(kron(kron(H_z, h_y*speye(m_y)), h_x*speye(m_x))))/3;
        H = sparse(kron(kron(H_y,H_x),H_z));
        
        
        % H-norm
        H_norm(i) = sqrt(diff'*H*diff);
%         H_norm(i) = sqrt(diff'*H*diff*h(i)^3);
        
        % L2-norm
        l2_norm(i) = sqrt(sum(diff.^2)*h(i)^3);
        
    end

    disp(['L2-norms: ', num2str(l2_norm)]);
    disp(['H-norms: ', num2str(H_norm)]);
    
    figure('Name', 'Convergence study');
    loglog(h, l2_norm, 'or');
    hold on
    loglog(h, H_norm, 'ob');
    xlabel('h');
    ylabel('||e||^2');
    title('Logarithmic plot over the stepsize h and the L2- norm of the error e');
    k = l2_norm(5)/h(5)^1;
    X = linspace(0.01, 0.2, 100);
    Y = k.*X.^1;
    loglog(X, Y, '-.k');
    k = H_norm(5)/h(5)^1;
    X = linspace(0.01, 0.2, 100);
    Y = k.*X.^1;
    loglog(X, Y, '-.k');
%     for i = 1:4
%         k = H_norm(5)/h(5)^i;
%         X = linspace(0.01, 0.2, 100);
%         Y = k.*X.^i;
%         loglog(X, Y, '-.k');
%     end
end