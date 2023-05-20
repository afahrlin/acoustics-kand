
hs = 0.1:0.01:0.4;
% max_es = zeros(1, length(hs));
% max_es_R = zeros(1, length(hs));
% 
% for i = 1:length(hs)
%     h = hs(i);
%     m_x = round(6.24/h);
%     m_y = round(3.84/h);
%     m_z = round(2.4/h);
%     dim = [m_x m_y m_z];
% 
%     R = 0;
%     max_es(i) = simulation_3D_stability(dim, R);
%     disp(['h = ', num2str(hs(i)), ' Done 1/2'])
%     R = 0.7;
%     max_es_R(i) = simulation_3D_stability(dim, R);
%     disp(['h = ', num2str(hs(i)), ' Done'])
% end

load('max_es')
load('max_es_R')

% Fit quadratic line to the data
coeffs = polyfit(hs, 1./max_es_R, 2);  % Use 2 for quadratic fit
disp(coeffs)

% Generate fitted values
h_fit = linspace(0*min(hs), max(hs), 100);  % Generate evenly spaced values of 'a' for fitting
es_fit = polyval(coeffs, h_fit);  % Evaluate the fitted quadratic line at the values of 'a'

disp(mean(max_es./max_es_R))

% figure(1)
% 
% plot(hs, 1./max_es, 'ko')
% hold on
% title('Maximum Eigenvalue at Different Grid Step Sizes h')
% plot(h_fit,f(h_fit),'k');
% %plot(h_fit,es_fit,'k');
% 
% plot(hs, 1./max_es_R, 'k+')
% plot(h_fit,f2(h_fit),'k');
% %plot(h_fit,1.85*es_fit,'k');
% xlabel('Grid Step Size h')
% xlim([0*min(hs) max(hs)])
% ylabel('Inverted Maximum Magnitude Eigenvalue 1/|\lambda|')
% % legend('Recorded Maximum Magnitude Eigenvalue', ...
% %     ['Fitted function: \lambda = 1/(', num2str(round(coeffs(1)*10000,2)), ...
% %     'e-4*h-', num2str(round(coeffs(2)*1000000,2)), 'e-6)'])
% % legend('Calculated Maximum Magnitude Eigenvalue, R = 0', ...
% %     'Fitted function: \lambda = 1/(4.5e-4*h - 4.5e-6)', ...
% %     'Calculated Maximum Magnitude Eigenvalue, R = 0.7', ...
% %     'Fitted function: \lambda = 1/(8.3e-4*h - 8.3e-6)')
% legend('Calculated Maximum Magnitude Eigenvalue, R = 0', ...
%     'Fitted function: 1/\lambda = 10^{-4}(1.556h^2+3.866h)', ...
%     'Calculated Maximum Magnitude Eigenvalue, R = 0.7', ...
%     'Fitted function: 1/\lambda = 10^{-4}(2.879h^2+7.152h)')
% 
% % save('max_es')
% % save('max_es_R')
% 
% function y = f(x)
% %     y = 1e-3*(0.155647852041533*x.^2+0.386592305792715*x+0.001259419760354);
%     y = 1e-3*(0.155647852041533*x.^2+0.386592305792715*x);
% 
% %     y = 4.5e-4*x - 4.5e-6;
% %     y = 4.2e-4*x;
% %     y = 1./(4.5e-4*x - 4.5e-6);
% %     y = 1./((4.5e-4)*x);
% end
% 
% function y = f2(x)
% %     y = 1/8.3e-4*1./(x - 0.01);
%     y = f(x)*1.85;
% %     y = 4.5e-4*x - 4.5e-6;
% %     y = 1e-3*x - 0.01e-3;
% end


