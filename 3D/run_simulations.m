

freqs = linspace(195, 300, 22);

for i = freqs
    disp(['Running frequency: ', num2str(i)]);
    simulation_3D_fourth(i);
end

disp(['Frequencies ', num2str(freqs(1)), ' to ', num2str(freqs(11)), ' done.'])


% sims = ["30Hz_9929", '35Hz_9582', '40Hz_4989', '45Hz_6189', '50Hz_7774', ...
%     '55Hz_6271', '60Hz_3118', '65Hz_7391', '70Hz_4478', '75Hz_2556',  ...
%     '80Hz_7737', '85Hz_6225', '90Hz_9463', '95Hz_7357', '100Hz_9952', ...
%     '105Hz_2383', '110Hz_8394', '115Hz_2365', '120Hz_4865', '125Hz_9377', ...
%     '130Hz_4611', '135Hz_5892', '140Hz_6514', '145Hz_2835', '150Hz_2636', ...
%     '155Hz_7751', '160Hz_3928', '165Hz_5914', '170Hz_1917', '175Hz_8814', ...
%     '180Hz_3849', '185Hz_2322', '190Hz_8652'];
% 
% disp(sims)
% 
% for i = sims
%     disp(['Running frequency: ', num2str(i)]);
%     plot_average(i);
% end
% 
% disp(['Frequencies ', num2str(sims(1)), ' to ', num2str(sims(end)), ' done.'])