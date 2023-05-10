% ========================================
% Operators
% 
% Based on code from the course scientific 
% computing for partial differential 
% equations, uppsala university
% ========================================

function [H, HI, D1, D2, e_l, e_r, d1_l, d1_r] = sbp_cent_4th(m, h)
    e_l = zeros(1, m);
    e_l(1) = 1;

    e_r = zeros(1, m);
    e_r(end) = 1;

    v = [17/48,59/48,43/48,49/48];

    H = eye(m);
    H(1:4, 1:4) = (diag(v));
    H(end-3:end, end-3:end) = (diag(flip(v)));
    H = H*h;

    HI = inv(H);
    
    %x1 = 0.70127127127127;
    
    Q = -1/12*diag(ones(m-2,1),2) + 8/12*diag(ones(m-1,1),1) - 8/12*diag(ones(m-1,1),-1) + 1/12*diag(ones(m-2,1),-2);
    Q_U = [0, 59/96, -1/12, -1/32; 
        -59/96, 0, 59/96, 0;
        1/12, -59/96, 0, 59/96;
        1/32, 0, -59/96, 0];
    Q(1:4,1:4) = Q_U;
    Q(end-3:end,end-3:end) = flip(flip(-Q_U,1),2);
    D1 = HI*(Q - 0.5*e_l'*e_l + 0.5*e_r'*e_r);
    
    M = 1/12*diag(ones(m-2,1),2) - 16/12*diag(ones(m-1,1),1) - 16/12*diag(ones(m-1,1),-1) + 1/12*diag(ones(m-2,1),-2) + 30/12*diag(ones(m,1),0);
    M_U = [9/8, -59/48, 1/12, 1/48; 
        -59/48, 59/24, -59/48, 0;
        1/12, -59/48, 55/24, -59/48;
        1/48, 0, -59/48, 59/24];
    
    M(1:4, 1:4) = M_U;
    M(end-3:end, end-3:end) = flip(flip(M_U,1),2);
    M = M/h;

    d_stenc = [-11/6, 3, -3/2, 1/3]/h;
    d1_l = zeros(1, m);
    d1_l(1:4) = d_stenc;
    d1_r = zeros(1, m);
    d1_r(end-3:end) = flip(-d_stenc);
    
    D2 = HI*(-M - e_l'*d1_l + e_r'*d1_r);
    disp(D2);
end







