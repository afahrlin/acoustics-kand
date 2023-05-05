% Convergence study main

function main()
    [uref, id] = make_ref(100);
    location = append('Testdata/', id, '/');
    save(location, 'uref', '-v7.3');
end

