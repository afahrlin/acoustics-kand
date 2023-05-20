% Code to plot stability region of rk4 method

x = linspace(-3, 1, 100);
y = linspace(-3, 3, 100);
[X, Y] = meshgrid(x,y);
Z = X + i*Y;
R = stability(Z);

x2 = linspace(-3, 0, 100);
y2 = linspace(-3, 3, 100);
[X2, Y2] = meshgrid(x2,y2);
Z2 = X2 + i*Y2;
C = f(Z2);
contour(x,y,abs(R),[1 1],'k')
hold on
contour(x2,y2,abs(C),[1 1],'--k')
title('Absolute Stability of RK4')
xlabel('Re(\lambdah_t)')
ylabel('Im(\lambdah_t)')
xline(0)
yline(0)
pbaspect([4 6 1])


function R = stability(z)
    R = 1 + z + 1/2*z.^2 + 1/6*z.^3 + 1/24*z.^4;
end

function res = f(z)
    res = abs(z/2.6);
end


