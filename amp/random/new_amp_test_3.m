% Make data to find actual value of coeficcient in distance formula
function new_amp_test_3()
    dim = [573,353,221];
    f = 200;
    w = 2*pi*f;
    T = 1/343 + pi/(2*w);

    disp('Testing with alva = 1/17.0852');
    disp('Dimensions: 11');
    [u1, simname, alva] = make_ref(f, T, dim, 1/17.0852, '1by17-11');
    location = append('Testdata/', simname, '/');

    save(append(location, 'u1.mat'), 'u1', '-v7.3');
    disp('Data saved');
    disp(simname);
    disp(['Alva: ', alva]);
end