
function memory_test()

all_points = sortrows(cell2mat(struct2cell(load('points.mat'))));

f = 200;
T = 0.001;
olle = 1;
mkdir('Testdata');

for i = 3:length(all_points)
    dim = all_points(i,:);
    [u, simname, alva] = make_ref(f, T, dim, olle);
    disp(alva);
    disp(' ');
end

end