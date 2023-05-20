
load('amp/INFO.mat');
u = cell2mat(struct2cell(load('amp/u1.mat', 'u1')));
load('amp/INFO.mat');
size(u)
%alva = u(101, 193, 121) % refgrid
alva = [u(91, 177, 111), u(92, 177, 111), u(93, 177, 111)]; % 11-grid

