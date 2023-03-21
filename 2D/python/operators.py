import numpy as np
import scipy.linalg as splg
import scipy.sparse as spsp

# Central 1D second order accurate finite difference SBP operators.
# Input:
#   m - number of grid points (integer)
#   h - step size (float)
# 
# Output:
#   H - inner product matrix
#   HI - inverse of H
#   D1 - first derivative SBP operator
#   D2 - second derivative SBP operator
#   e_l,e_r - vectors to extract the boundary grid points
#   d1_l,d1_r - vectors to extract the first derivatives at the boundary grid points
# 
# Use as follows:
# 
# import operators as ops
# H,HI,D1,D2,e_l,e_r,d1_l,d1_r = ops.sbp_cent_2nd(m,h,order)
def sbp_cent_2nd(m,h):
    e_l = np.zeros(m)
    e_l[0] = 1

    e_r = np.zeros(m)
    e_r[-1] = 1

    H = np.eye(m)
    H[0,0] = 0.5
    H[-1,-1] = 0.5
    H = h*H

    HI = np.linalg.inv(H)

    D1 = 0.5*np.diag(np.ones(m-1),1) - 0.5*np.diag(np.ones(m-1),-1)
    D1[0,0] = -1
    D1[0,1] = 1
    D1[-1,-2] = -1
    D1[-1,-1] = 1
    D1 = D1/h

    Q = np.matmul(H,D1) + 0.5*np.tensordot(e_l, e_l, axes=0) - 0.5*np.tensordot(e_r, e_r, axes=0)

    D2 = np.diag(np.ones(m-1),1) + np.diag(np.ones(m-1),-1) - 2*np.diag(np.ones(m),0)
    D2[0,0] = 1
    D2[0,1] = -2
    D2[0,2] = 1
    D2[-1,-3] = 1
    D2[-1,-2] = -2
    D2[-1,-1] = 1
    D2 = D2/(h*h)

    d_stenc = np.array([-3./2, 2, -1./2])/h
    d1_l = np.zeros(m)
    d1_l[:3] = d_stenc
    d1_r = np.zeros(m)
    d1_r[-3:] = -np.flip(d_stenc)

    M = -np.matmul(H,D2) - np.tensordot(e_l, d1_l, axes=0) + np.tensordot(e_r, d1_r, axes=0)

    H = spsp.csc_matrix(H)
    HI = spsp.csc_matrix(HI)
    D1 = spsp.csc_matrix(D1)
    D2 = spsp.csc_matrix(D2)
    e_l = spsp.csc_matrix(e_l)
    e_r = spsp.csc_matrix(e_r)
    d1_l = spsp.csc_matrix(d1_l)
    d1_r = spsp.csc_matrix(d1_r)
    return H,HI,D1,D2,e_l,e_r,d1_l,d1_r

# Central 1D fourth order accurate finite difference SBP operators.
# Input:
#   m - number of grid points (integer)
#   h - step size (float)
# 
# Output:
#   H - inner product matrix
#   HI - inverse of H
#   D1 - first derivative SBP operator
#   D2 - second derivative SBP operator
#   e_l,e_r - vectors to extract the boundary grid points
#   d1_l,d1_r - vectors to extract the first derivatives at the boundary grid points
# 
# Use as follows:
# 
# import operators as ops
# H,HI,D1,D2,e_l,e_r,d1_l,d1_r = ops.sbp_cent_4th(m,h,order)
def sbp_cent_4th(m,h):
    e_l = np.zeros(m)
    e_l[0] = 1

    e_r = np.zeros(m)
    e_r[-1] = 1

    H = np.diag(np.ones(m))
    H[0:4,0:4] = np.diag(np.array([17/48, 59/48, 43/48, 49/48]))
    H[-4:,-4:] = np.diag(np.array([49/48, 43/48, 59/48, 17/48]))
    H=H*h;

    HI = np.linalg.inv(H)

    Q = -1/12*np.diag(np.ones(m-2),2) + 8/12*np.diag(np.ones(m-1),1) - 8/12*np.diag(np.ones(m-1),-1) + 1/12*np.diag(np.ones(m-2),-2)
    Q_U = np.array([[0, 0.59e2/0.96e2, -0.1e1/0.12e2, -0.1e1/0.32e2],[-0.59e2/0.96e2, 0, 0.59e2/0.96e2, 0],[0.1e1/0.12e2, -0.59e2/0.96e2, 0, 0.59e2/0.96e2],[0.1e1/0.32e2, 0, -0.59e2/0.96e2, 0]])
    Q[0:4,0:4] = Q_U;
    Q[-4:,-4:] = np.flipud(np.fliplr(-Q_U));

    D1 = HI@(Q - 0.5*np.tensordot(e_l, e_l, axes=0) + 1/2*np.tensordot(e_r, e_r, axes=0))


    M_U = np.array([[0.9e1/0.8e1, -0.59e2/0.48e2, 0.1e1/0.12e2, 0.1e1/0.48e2],[-0.59e2/0.48e2, 0.59e2/0.24e2, -0.59e2/0.48e2, 0],[0.1e1/0.12e2, -0.59e2/0.48e2, 0.55e2/0.24e2, -0.59e2/0.48e2],[0.1e1/0.48e2, 0, -0.59e2/0.48e2, 0.59e2/0.24e2]])
    M = -(-1/12*np.diag(np.ones(m-2),2) + 16/12*np.diag(np.ones(m-1),1) + 16/12*np.diag(np.ones(m-1),-1) - 1/12*np.diag(np.ones(m-2),-2) - 30/12*np.diag(np.ones(m),0));

    M[0:4,0:4] = M_U

    M[-4:,-4:] = np.flipud(np.fliplr(M_U))
    M=M/h;

    d_stenc = np.array([-0.11e2/0.6e1, 3, -0.3e1/0.2e1, 0.1e1/0.3e1])/h
    d1_l = np.zeros(m)
    d1_l[0:4] = d_stenc
    d1_r = np.zeros(m)
    d1_r[-4:] = np.flip(-d_stenc)

    D2 = HI@(-M - np.tensordot(e_l, d1_l, axes=0) + np.tensordot(e_r, d1_r, axes=0))

    H = spsp.csc_matrix(H)
    HI = spsp.csc_matrix(HI)
    D1 = spsp.csc_matrix(D1)
    D2 = spsp.csc_matrix(D2)
    e_l = spsp.csc_matrix(e_l)
    e_r = spsp.csc_matrix(e_r)
    d1_l = spsp.csc_matrix(d1_l)
    d1_r = spsp.csc_matrix(d1_r)
    return H,HI,D1,D2,e_l,e_r,d1_l,d1_r

# Central 1D sixth order accurate finite difference SBP operators.
# Input:
#   m - number of grid points (integer)
#   h - step size (float)
# 
# Output:
#   H - inner product matrix
#   HI - inverse of H
#   D1 - first derivative SBP operator
#   D2 - second derivative SBP operator
#   e_l,e_r - vectors to extract the boundary grid points
#   d1_l,d1_r - vectors to extract the first derivatives at the boundary grid points
# 
# Use as follows:
# 
# import operators as ops
# H,HI,D1,D2,e_l,e_r,d1_l,d1_r = ops.sbp_cent_6th(m,h,order)
def sbp_cent_6th(m,h):
    e_l = np.zeros(m)
    e_l[0] = 1

    e_r = np.zeros(m)
    e_r[-1] = 1

    H = np.diag(np.ones(m),0);
    H[:6,:6] = np.diag(np.array([13649/43200,12013/8640,2711/4320,5359/4320,7877/8640, 43801/43200]))
    H[-6:,-6:] = np.fliplr(np.flipud(np.diag(np.array([13649/43200,12013/8640,2711/4320,5359/4320,7877/8640,43801/43200]))));
    H=H*h;

    HI = np.linalg.inv(H)

    x1 = 0.70127127127127;

    D1 = 1/60*np.diag(np.ones(m-3),3) - 9/60*np.diag(np.ones(m-2),2) + 45/60*np.diag(np.ones(m-1),1) - 45/60*np.diag(np.ones(m-1),-1) + 9/60*np.diag(np.ones(m-2),-2) - 1/60*np.diag(np.ones(m-3),-3)

    D1_bound_stencil = np.array([[-21600/13649, 43200/13649*x1-7624/40947, -172800/13649*x1 + 715489/81894, 259200/13649*x1-187917/13649, -172800/13649*x1+735635/81894, 43200/13649*x1-89387/40947, 0, 0, 0], \
        [-8640/12013*x1+7624/180195, 0, 86400/12013*x1-57139/12013, -172800/12013*x1+745733/72078, 129600/12013*x1-91715/12013,-34560/12013*x1+240569/120130, 0, 0, 0], \
        [17280/2711*x1-715489/162660, -43200/2711*x1+57139/5422, 0, 86400/2711*x1-176839/8133, -86400/2711*x1+242111/10844, 25920/2711*x1-182261/27110, 0, 0, 0], \
        [-25920/5359*x1+187917/53590, 86400/5359*x1-745733/64308, -86400/5359*x1+176839/16077, 0, 43200/5359*x1-165041/32154, -17280/5359*x1+710473/321540, 72/5359, 0, 0], \
        [ 34560/7877*x1-147127/47262, -129600/7877*x1+91715/7877, 172800/7877*x1-242111/15754, -86400/7877*x1+165041/23631, 0, 8640/7877*x1, -1296/7877, 144/7877, 0], \
        [-43200/43801*x1+89387/131403, 172800/43801*x1-240569/87602, -259200/43801*x1+182261/43801, 172800/43801*x1-710473/262806, -43200/43801*x1, 0, 32400/43801, -6480/43801, 720/43801]])

    D1[:6,:9] = D1_bound_stencil
    D1[-6:,-9:] = np.flipud(np.fliplr(-D1_bound_stencil))
    D1 = D1/h

    Q = np.matmul(H,D1) + 0.5*np.tensordot(e_l, e_l, axes=0) - 0.5*np.tensordot(e_r, e_r, axes=0)

    D2 = (2*np.diag(np.ones(m-3),3) - 27*np.diag(np.ones(m-2),2) + 270*np.diag(np.ones(m-1),1) + 270*np.diag(np.ones(m-1),-1) - 27*np.diag(np.ones(m-2),-2) + 2*np.diag(np.ones(m-3),-3) - 490*np.diag(np.ones(m),0))/180;

    D2_bound_stencil = np.array([[114170/40947, -438107/54596, 336409/40947, -276997/81894, 3747/13649, 21035/163788, 0, 0, 0], \
        [6173/5860, -2066/879, 3283/1758, -303/293, 2111/3516, -601/4395, 0, 0, 0], \
        [-52391/81330, 134603/32532, -21982/2711, 112915/16266, -46969/16266, 30409/54220, 0, 0, 0], \
        [68603/321540, -12423/10718, 112915/32154, -75934/16077, 53369/21436, -54899/160770, 48/5359, 0, 0], \
        [-7053/39385, 86551/94524, -46969/23631, 53369/15754, -87904/23631, 820271/472620, -1296/7877, 96/7877, 0], \
        [21035/525612, -24641/131403, 30409/87602, -54899/131403, 820271/525612, -117600/43801, 64800/43801, -6480/43801, 480/43801]]);

    D2[:6,:9] = D2_bound_stencil
    D2[-6:,-9:] = np.flipud(np.fliplr(D2_bound_stencil))

    D2 = D2/h**2

    d_stenc = np.array([-25/12, 4, -3, 4/3, -1/4])/h
    d1_l = np.zeros(m)
    d1_l[0:5] = d_stenc
    d1_r = np.zeros(m)
    d1_r[-5:] = np.flip(-d_stenc)

    M = -np.matmul(H,D2) - np.tensordot(e_l, d1_l, axes=0) + np.tensordot(e_r, d1_r, axes=0)

    H = spsp.csc_matrix(H)
    HI = spsp.csc_matrix(HI)
    D1 = spsp.csc_matrix(D1)
    D2 = spsp.csc_matrix(D2)
    e_l = spsp.csc_matrix(e_l)
    e_r = spsp.csc_matrix(e_r)
    d1_l = spsp.csc_matrix(d1_l)
    d1_r = spsp.csc_matrix(d1_r)
    return H,HI,D1,D2,e_l,e_r,d1_l,d1_r

