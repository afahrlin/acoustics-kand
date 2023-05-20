% test convergence based on UPPMAX runs
% H-norm

% Load everything
function conv_study()
    uref = struct2array(load('Konvergens/REF/u.mat', 'u'));
    % normalise uref (illegal?)
%     uref = uref./max(max(max(uref)));

%     u1 = struct2array(load(append('Testdata/conv_test/200Hz_36729points__1/u.mat'), 'u'));
%     u2 = struct2array(load(append('Testdata/conv_test/200Hz_279825points__2/u.mat'), 'u'));
%     u3 = struct2array(load(append('Testdata/conv_test/200Hz_928969points__3/u.mat'), 'u'));
%     u4 = struct2array(load(append('Testdata/conv_test/200Hz_2183841points__4/u.mat'), 'u'));
%     u6 = struct2array(load(append('Testdata/conv_test/200Hz_7309489points__6/u.mat'), 'u'));
%     
    u1 = struct2array(load(append('Konvergens/t1/u.mat'), 'u'));
    u2 = struct2array(load(append('Konvergens/t2/u.mat'), 'u'));
    u3 = struct2array(load(append('Konvergens/t3/u.mat'), 'u'));
    u4 = struct2array(load(append('Konvergens/t4/u.mat'), 'u'));
    u6 = struct2array(load(append('Konvergens/t6/u.mat'), 'u'));
   
    n = [12, 6, 4, 3, 2];
    h = [0.12, 0.06, 0.04, 0.03, 0.02];
    all_norm = zeros(1, 5);

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
        
        % normalise u (illegal?)
        %u = u./max(max(max(u)));
        
        % parameters relevant to this size
        ref = uref(1:n(i):end, 1:n(i):end, 1:n(i):end);
        m_x = size(u, 1);
        m_y = size(u, 2);
        m_z = size(u, 3);
        
        % make H
        v = [17/48,59/48,43/48,49/48];
        
        H_x = speye(m_x);
        H_x(1:4, 1:4) = (diag(v));
        H_x(end-3:end, end-3:end) = (diag(flip(v)));
        H_x = H_x*h(i);
        
        H_y = speye(m_y);
        H_y(1:4, 1:4) = (diag(v));
        H_y(end-3:end, end-3:end) = (diag(flip(v)));
        H_y = H_y*h(i);
        
        H_z = speye(m_z);
        H_z(1:4, 1:4) = (diag(v));
        H_z(end-3:end, end-3:end) = (diag(flip(v)));
        H_z = H_z*h(i);
        
        H_xy = sparse(kron(speye(m_y), H_x) + kron(H_y, speye(m_x)));
        H = sparse(kron(speye(m_z), H_xy) + kron(H_z, speye(m_x*m_y)));
        clear H_xy H_z
        
        % difference
        diff = ref - u;
        diff = reshape(diff, numel(diff), 1);
        norm = diff' * H * diff;

        all_norm(i) = norm;
        disp(norm);
        disp(i);
    end

    disp(all_norm);
end

