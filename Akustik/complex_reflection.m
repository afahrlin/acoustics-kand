
f = 20;

rho_a = 1.293e-3;
c_a = 343;

rho_w = 700;
c_w = 4570;

L_w = 0.05;
l_w = c_w/f;
k_w = 1/l_w;

R_aw = abs(reflection(rho_a, c_a, rho_w, c_w, L_w, f));
R_wa = abs(reflection(rho_w, c_w, rho_a, c_a, L_w, f));
T_aw = abs(transmission(rho_a, c_a, rho_w, c_w, L_w, f));
T_wa = abs(transmission(rho_w, c_w, rho_a, c_a, L_w, f));

R = abs(R_aw + T_aw*T_wa*R_wa*exp(2*1i*L_w*k_w) + T_aw*T_wa/(1-R_wa^2*exp(2*1i*L_w*k_w)))
R = abs(reflection(rho_a, c_a, rho_w, c_w, L_w, f))

function T = transmission(rho_a, c_a, rho_w, c_w, L_w, f)
    l_a = c_a/f;
    k_a = 1/l_a;

    l_w = c_w/f;
    k_w = 1/l_w;

    r_a = rho_a*c_a;
    r_w = rho_w*c_w;

    sqr_pos = (r_a+r_w)^2;
    sqr_neg = (r_a-r_w)^2;
    
    T = (4*(r_a*r_w/sqr_pos)*exp(1i*L_w*(k_w-k_a)))/(1-exp(2*1i*L_w*k_w)*sqr_neg/sqr_pos);
end

function R = reflection(rho_a, c_a, rho_w, c_w, L_w, f)
    l_w = c_w/f;
    k_w = 1/l_w;

    r_a = rho_a*c_a;
    r_w = rho_w*c_w;

    sqr_neg = (r_a-r_w)^2;
    sqr_pos = (r_a+r_w)^2;
    compl = exp(2*1i*L_w*k_w);

    R = (((r_a^2-r_w^2)*compl-r_a^2+r_w^2)/sqr_pos)/(1-compl*sqr_neg/sqr_pos);
end




