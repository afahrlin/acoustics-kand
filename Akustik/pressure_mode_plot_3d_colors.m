% Code to generate plots of mode (n_x,n_y)

L_x = 6.24;
L_y = 3.84;
L_z = 2.4;
n_x = 2;
n_y = 2;
n_z = 0;

x = linspace(0, L_x, 100);
y = linspace(0, L_y, 100);
z = linspace(0, L_z, 100);

[X, Y, Z] = meshgrid(x,y,z);
[X_plot_top, Y_plot_top] = meshgrid(x,y);
[X_plot_side, Z_plot_side] = meshgrid(x,z);
[Y_plot_front, Z_plot_front] = meshgrid(y,z);

c = [-1 1];

P = pressure_mode(1.00001, X, Y, Z, n_x, n_y, n_z, L_x, L_y, L_z);

% newpoints = 100;
% [xq,yq] = meshgrid(...
%             linspace(min(min(xrow,[],2)),max(max(xrow,[],2)),newpoints ),...
%             linspace(min(min(ycol,[],1)),max(max(ycol,[],1)),newpoints )...
%           );
% BDmatrixq = interp2(xrow,ycol,BDmatrix,xq,yq,'cubic');
% [c,h]=contourf(xq,yq,BDmatrixq);

P = permute(P, [2 1 3]);

figure();
v = -1:0.01:1;
v2 = [-1,-.75,-.5,-.25,0,.25,.50,.75,1];
[C,h] = contourf(X_plot_top, Y_plot_top, transpose(P(:,:,end)), v);
set(h, 'edgecolor','none');
hold on
[C,h] = contour(X_plot_top, Y_plot_top, transpose(P(:,:,end)), v2, 'k');
clabel(C,h,v2,'FontSize',8,'labelspacing', 1000)
% colormap jet
caxis(c);
%cb1 = colorbar;
%cb1.Label.String = 'Normalized Sound Pressure';
title(['Top view of pressure distribution of mode (', ...
    num2str(n_x), ',', num2str(n_y), ',', num2str(n_z), ')']);
xlabel('x');
ylabel('y');
pbaspect([L_x L_y, 1]);

figure();
[C,h] = contourf(X_plot_side, Z_plot_side, transpose(reshape(P(:,1,:), numel(x), numel(z))), v);
%colormap jet
set(h, 'edgecolor','none');
hold on
[C,h] = contour(X_plot_side, Z_plot_side, transpose(reshape(P(:,1,:), numel(x), numel(z))), v2, 'k');
clabel(C,h,v2,'FontSize',8,'labelspacing', 1000)
caxis(c);
cb2 = colorbar;
cb2.Label.String = 'Normalized Absolute Sound Pressure';
title(['Side view of pressure distribution of mode (', ...
    num2str(n_x), ',', num2str(n_y), ',', num2str(n_z), ')']);
xlabel('x');
ylabel('z');
pbaspect([L_x L_z, 1]);

% figure();
% [C,h] = contourf(Y_plot_front, Z_plot_front, transpose(reshape(P(1,:,:), numel(y), numel(z))), v);
% set(h, 'edgecolor','none');
% %colormap jet
% caxis(c);
% cb3 = colorbar;
% cb3.Label.String = 'Normalized Sound Pressure';
% title(['Front view of pressure distribution of mode (', ...
%     num2str(n_x), ',', num2str(n_y), ',', num2str(n_z), ')']);
% xlabel('y');
% ylabel('z');
% pbaspect([L_y L_z, 1]);

function p = pressure_mode(C, x, y, z, n_x, n_y, n_z, L_x, L_y, L_z)
    p = C*cos(n_x*pi*x/L_x).*cos(n_y*pi*y/L_y).*cos(n_z*pi*z/L_z);
end

