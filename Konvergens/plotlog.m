    
% plot pls

%y = [0.7173, 0.7496, 0.8407, 0.9042]*10^5;
% y = [4.0703, 3.2278, 2.5695, 2.0123, 1.1224];
% y = [208.7216  168.4450  167.1107  173.6287  176.9409];
% y = [1.2077    1.0593    0.9400    0.8295    0.6174];
% y = [25.0466   10.1067    6.6844    5.2089    3.5388]; % l2 not
% normalized
y = [0.1449    0.0636    0.0376    0.0249    0.0123]; % l2 normalised 
x = [0.12, 0.06, 0.04, 0.03, 0.02]; 

figure('Name', 'Convergence study');
loglog(x, y, 'ok');
xlabel('h');
ylabel('||e||^2');
title('Logarithmic plot over the stepsize h and the L2- norm of the error e');
hold on
k = y(5)/x(5)^2;
X = linspace(0.01, 0.2, 100);
Y = k.*X.^2;
loglog(X, Y, '-.k');
hold off
