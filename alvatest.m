
load('INFO.mat');
u = cell2mat(struct2cell(load('u1.mat', 'u1')));
U = reshape(u, numel(u), 1);
alva = U(round(m*0.5) - m_x*round(0.5*m_y) + m_y*round(0.5*m_x) + round(m_x/L_x));
disp(alva);
disp(1/alva);

