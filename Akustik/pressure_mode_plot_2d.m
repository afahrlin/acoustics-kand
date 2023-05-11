% Code to generate plots of mode (n_x,n_y)

L_x = 5;
L_y = 3;
n_x = 3;
n_y = 2;

x = linspace(0, L_x, 100);
y = linspace(0, L_y, 100);

[X, Y] = meshgrid(x,y);

P = pressure_mode(1.000000001, X, Y, n_x, n_y, L_x, L_y);

figure();
v = -1:0.25:1;
v2 = [-1,-.75,-.5,-.25,.25,.50,.75,1];
[C,h] = contour(X, Y, P, v, 'k');
clabel(C,h,v,'FontSize',8,'labelspacing', 1000, 'BackgroundColor', 'w')
title(['Pressure distribution of mode (', num2str(n_x), ',', num2str(n_y), ')']);
xlabel('x');
ylabel('y');

function p = pressure_mode(C, x, y, n_x, n_y, L_x, L_y)
    p = C*cos(n_x*pi*x/L_x).*cos(n_y*pi*y/L_y);
end

