% Convergence study main

function main()
    T = 0.3;
    f = 160;
    
    [uref, simname] = make_ref(f, T);
    location = append('Testdata/', simname);
    save(append(location, '/info2.mat'), 'f', 'T', '-v7.3');
    
    save(append(location, '/uref.mat'), 'uref', '-v7.3');
    disp('Data saved');
    disp(simname);
end
