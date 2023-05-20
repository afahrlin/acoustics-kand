% Make data to find actual value of coeficcient in distance formula
function new_amp_test()
    dim = [625, 385, 241];
    f = 200;
    w = 2*pi*f;
    T = 1/343 + pi/(2*w);

    [u1, simname, alva] = make_ref(f, T, dim, 1, 'original');
    location = append('Testdata/', simname, '/');

    save(append(location, 'u1.mat'), 'u1', '-v7.3');
    disp('Data saved');
    disp(simname);
    disp(['Alva: ', alva]);
end