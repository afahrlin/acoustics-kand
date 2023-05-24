% test convergence based on UPPMAX runs
% L2 norm

% Load everything
function conv_study2()
%     uref = struct2array(load('Konvergens/REF/u.mat', 'u'));
    uref = struct2array(load('../Testdata/conv_test/REF/u.mat', 'u'));
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

        % normalising u, (illegal?)
%         u = u./max(max(max(u)));

        % values relevant to this size
        ref = uref(1:n(i):end, 1:n(i):end, 1:n(i):end);

        % difference
        diff = ref - u;
        diff = reshape(diff, numel(diff), 1);
        diff2 = diff.*diff;
        norm = sqrt(sum(diff2)*h(i)^3);

        l2_norm(i) = norm;
        disp(norm);
        disp(i);
    end

    disp(l2_norm);
end