% ================
% Rungekutta 4
% ================

function [v, t] = rungekutta4(f, A, v, t, dt)
    k1 = dt*f(v, A);
    k2 = dt*f(v + 0.5*k1, A);
    k3 = dt*f(v + 0.5*k2, A);
    k4 = dt*f(v + k3, A);

    v = v + 1/6*(k1 + 2*k2 + 2*k3 + k4);
    t = t + dt;
end