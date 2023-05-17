% Convergence test for 3D wave equation

f = 70; % Hz
T = 0.4; % s

all_points = sortrows(cell2mat(struct2cell(load('Konvergens/points.mat'))));

location = 'Testdata/conv_test';
mkdir(location);

[uref, refname, ~] = make_ref(f, T, all_points(12));
save(append(location, '/ref.mat'), 'uref', 'refname', '-v7.3');
disp('Reference made');

for i = [1, 2, 3, 6]
    dim = all_points(i,1:3);
    n = 12/all_points(i, 4);
    
    [u, simname, H] = make_test(f, T, dim);
    
    simname = append(simname, '_', num2str(n));
    disp(append('Test made ', num2str(n))); 
    
    ref = permute(uref, [2,1,3]);        
    ref = ref(1:n:end, 1:n:end, 1:n:end);
    ref = reshape(ref, numel(ref), 1);
    
    diff = ref - reshape(u, numel(u), 1);
    norm = diff'*H*diff;
    disp(append('Norm: ', norm));
    save(append(location, '/', simname, '.mat'), 'u', 'diff', 'norm', 'n', 'H', '-v7.3');
    
end

disp('All tests made');

for i = [1, 2, 3, 6]
    n = all_points(i, 4);
    
    ref = permute(uref, [2,1,3]);        
    ref = ref(1:n:end, 1:n:end, 1:n:end);
    ref = reshape(ref, numel(ref), 1);
    
    diff = ref - 

end

