% Convergence test for 3D wave equation

f = 70; % Hz
T = 0.4; % s

load('Konvergens/points.mat');

[uref, refname] = make_ref(f, T, k5);
location = append('Testdata/', refname);
save(append(location, '/conv_test.mat'), 'uref', '-v7.3');

[u1, name_1] = make_test(f, T, k1);
[u2, name_2] = make_test(f, T, k2);
[u3, name_3] = make_test(f, T, k3);
[u4, name_4] = make_test(f, T, k4);

save(append(location, '/conv_test.mat'), 'u1', 'u2', 'u3', 'u4', '-v7.3');

disp('Data saved');
disp(simname);

% tests = [u1, u2, u3, u4];
% names = [name_1, name_2, name_3, name_4];
% max_diffs = zeros(1, 4);
%
% for i = 1:4
%     s = 2^(5 - i);
%     newref = uref(1:s:end, 1:s:end, 1:s:end);
%     diff = tests(1) - newref;

%   MAX DIFF APPROACH
%     max_diff = max(abs(diff));
%     diffs(i) = max_diff;

%   L2 NORM APPROACH?
% end

