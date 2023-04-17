

n = 20;
amps = zeros(n, 1);
a = zeros(n, 1);

%for i = 0:1/(n-1):1
for i = 1:n+1
    a(i) = (i-1)/n;
    amps(i) = absorption_test(a(i));
end

% Fit quadratic line to the data
coeffs = polyfit(a, amps, 3);  % Use 2 for quadratic fit
disp(coeffs)

% Generate fitted values
a_fit = linspace(min(a), max(a), 100);  % Generate evenly spaced values of 'a' for fitting
amps_fit = polyval(coeffs, a_fit);  % Evaluate the fitted quadratic line at the values of 'a'

plot(a, amps, 'ko')
hold on
plot(a_fit,amps_fit,'k');
title('Absorption-rate of a One-Dimensional Wave')
xlabel('\beta_3')
ylabel('Amplitude of reflected wave')
ylim([0 1])
legend('Amplitude of reflected waves', 'Fitted function: y = -0.454x^3+1.340x^2-1.883x+0.994')


