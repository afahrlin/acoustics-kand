% Make data to find actual value of coeficcient in distance formula

dim = [833, 513, 321];
f = 200;
w = 2*pi*f;
T = 1/343 + pi/(2*w);

[uref, simname, alva] = make_ref(f, T, dim, 1);
location = append('Testdata/');
save(append(location, '/info2.mat'), 'f', 'T', '-v7.3');

save(append(location, '/uref.mat'), 'uref', '-v7.3');
disp('Data saved');
disp(simname);
disp(append('Alva: ', alva));

[u2, simname, alva] = make_ref(f, T, dim, 1/alva);

save(append(location, '/u2.mat'), 'u2', '-v7.3');
disp('Data saved');
disp(simname);
disp(append('Alva: ', alva));