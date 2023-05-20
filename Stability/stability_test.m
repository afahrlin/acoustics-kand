colors = ['r', 'g', 'c', 'b', 'm'];
hs = [0.2 0.225 0.25 0.275 0.3];
max_es = zeros(1, length(hs));

for c = 1:length(colors)
    h = hs(c);
    m_x = round(6.24/h);
    m_y = round(3.84/h);
    m_z = round(2.4/h);
    dim = [m_x m_y m_z];

    R = 0:0.05:1;
    es = zeros(1,length(R));

    for i = 1:length(R)
        es(i) = simulation_3D_stability(dim, R(i));
        disp(['R = ', num2str(R(i)), ' Done'])
    end
    
    max_es(c) = es(1);
    
    figure(1);
    plot(R, es, append(colors(c),'o'))
    hold on
    disp(['\n h = ', num2str(hs(c)), ' Done'])
end

title('Maximum Magnitude Eigenvalue at Different Reflection Rates')
ylim([0 inf])
xlabel('Reflection Rate R')
ylabel('Maximum Magnitude Eigenvalue')
legend(['h = ', num2str(hs(1))], ['h = ', num2str(hs(2))], ...
    ['h = ', num2str(hs(3))], ['h = ', num2str(hs(4))], ...
    ['h = ', num2str(hs(5))])




