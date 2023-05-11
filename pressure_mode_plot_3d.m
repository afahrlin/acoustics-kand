% Code to generate plots of mode (n_x,n_y)

L_x = 5;
L_y = 3;
L_z = 2;
n_x = 3;
n_y = 0;
n_z = 0;

x = linspace(0, L_x, 100);
y = linspace(0, L_y, 100);
z = linspace(0, L_z, 100);

[X, Y, Z] = meshgrid(x,y,z);
[X_plot_top, Y_plot_top] = meshgrid(x,y);
[X_plot_side, Z_plot_side] = meshgrid(x,z);
[Y_plot_front, Z_plot_front] = meshgrid(y,z);

P = pressure_mode(1.00001, X, Y, Z, n_x, n_y, n_z, L_x, L_y, L_z);

disp(size(X_plot_side))
disp(size(Z_plot_side))
disp(size(reshape(P(:,1,:), numel(x), numel(z))))

P = permute(P, [2 1 3]);

figure();
v = -1:0.25:1;
%v2 = [-1,-.75,-.5,-.25,.25,.50,.75,1];
[C,h] = contour(X_plot_top, Y_plot_top, transpose(P(:,:,1)), v, 'k');
clabel(C,h,v,'FontSize',8,'labelspacing', 1000, 'BackgroundColor', 'w')
title(['Top view of pressure distribution of mode (', ...
    num2str(n_x), ',', num2str(n_y), ',', num2str(n_z), ')']);
xlabel('x');
ylabel('y');
pbaspect([L_x L_y, 1]);

figure();
[C,h] = contour(X_plot_side, Z_plot_side, transpose(reshape(P(:,1,:), numel(x), numel(z))), v, 'k');
clabel(C,h,v,'FontSize',8,'labelspacing', 1000, 'BackgroundColor', 'w')
title(['Side view of pressure distribution of mode (', ...
    num2str(n_x), ',', num2str(n_y), ',', num2str(n_z), ')']);
xlabel('x');
ylabel('z');
pbaspect([L_x L_z, 1]);

figure();
[C,h] = contour(Y_plot_front, Z_plot_front, transpose(reshape(P(1,:,:), numel(y), numel(z))), v, 'k');
clabel(C,h,v,'FontSize',8,'labelspacing', 1000, 'BackgroundColor', 'w')
title(['Front view of pressure distribution of mode (', ...
    num2str(n_x), ',', num2str(n_y), ',', num2str(n_z), ')']);
xlabel('y');
ylabel('z');
pbaspect([L_y L_z, 1]);

function p = pressure_mode(C, x, y, z, n_x, n_y, n_z, L_x, L_y, L_z)
    p = C*cos(n_x*pi*x/L_x).*cos(n_y*pi*y/L_y).*cos(n_z*pi*z/L_z);
end

