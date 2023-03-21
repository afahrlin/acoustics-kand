import numpy as np
from scipy.sparse import kron, csc_matrix, eye, vstack
from math import ceil
import operators as ops
import matplotlib.pyplot as plt
from vispy import app, scene, gloo, visuals
import rungekutta4 as rk4
import time

import matplotlib
matplotlib.use('TkAgg')


def main():
    # Method parameters
    show_time_steps = False   # Set to False to only plot final values
    start_time = time.time()

    # Number of grid points
    m_x = 99
    m_y = 99
    m = (m_x+1)*(m_y+1)

    # Order of accuracy: 2, 4, or 6.
    order = 6

    # Model parameters
    T = 3           # Final time
    # T = 2/100
    # x_l = -3.43/2   # Left boundary of x
    # x_r = 3.43/2    # Right boundary of x
    x_l = -1#-1/2
    x_r = 1#1/2
    L_x = x_r-x_l   # Length of x interval
    # y_l = -3.43/2   # Left boundary of y
    # y_r = 3.43/2    # Right boundary of y
    y_l = -1/2
    y_r = 1/2
    L_y = y_r-y_l   # Length of y interval

    # PDE parameters
    # c = 343
    c = 1
    x_0 = 0
    y_0 = 0
    o = 0.05
    k_x = 2*np.pi
    k_y = 2*np.pi
    w = c*np.sqrt(k_x**2+k_y**2)
    
    # Spatial discretization
    h_x = L_x/(m_x)
    x_vec = np.linspace(x_l, x_r, m_x+1)
    h_y = L_y/(m_y)
    y_vec = np.linspace(y_l, y_r, m_y+1)
    X_vec, Y_vec = np.meshgrid(x_vec, y_vec)

    # Time discretization
    h_t = 0.1*h_x/c
    m_t = int(ceil(T/h_t))
    h_t = T/m_t
    
    # Get operators and default vectors from operators.py from lab1
    # Get D2 operator for x
    if order == 2:
        H_x, HI_x, D1_x, D2_x, e_l_x, e_r_x, d1_l_x, d1_r_x = ops.sbp_cent_2nd(m_x+1, h_x)
    elif order == 4:
        H_x, HI_x, D1_x, D2_x, e_l_x, e_r_x, d1_l_x, d1_r_x = ops.sbp_cent_4th(m_x+1, h_x)
    elif order == 6:
        H_x, HI_x, D1_x, D2_x, e_l_x, e_r_x, d1_l_x, d1_r_x = ops.sbp_cent_6th(m_x+1, h_x)
    else:
        raise NotImplementedError('Order not implemented.')
    
    # Get D2 operator for y
    if order == 2:
        H_y, HI_y, D1_y, D2_y, e_l_y, e_r_y, d1_l_y, d1_r_y = ops.sbp_cent_2nd(m_y+1, h_y)
    elif order == 4:
        H_y, HI_y, D1_y, D2_y, e_l_y, e_r_y, d1_l_y, d1_r_y = ops.sbp_cent_4th(m_y+1, h_y)
    elif order == 6:
        H_y, HI_y, D1_y, D2_y, e_l_y, e_r_y, d1_l_y, d1_r_y = ops.sbp_cent_6th(m_y+1, h_y)
    else:
        raise NotImplementedError('Order not implemented.')
    
    # SBP-SAT approximation, x
    D_x = c**2*D2_x + c**2*HI_x@e_l_x.transpose()@d1_l_x*eye(m_x+1) - c**2*HI_x@e_r_x.transpose()@d1_r_x*eye(m_x+1)
    # SBP-SAT approximation, y
    D_y = c**2*D2_y + c**2*HI_y@e_l_y.transpose()@d1_l_y*eye(m_y+1) - c**2*HI_y@e_r_y.transpose()@d1_r_y*eye(m_y+1)
    # SBP-SAT aprroximation
    D = (kron(eye(m_y+1), D_x)+kron(D_y, eye(m_x+1)))
    
    # Construct matrix: u_t = Au with u = [phi, phi_t]^T
    # [[0, I],
    #  [D, 0]]
    A = csc_matrix((2*m, 2*m))
    A[:m, m:] = eye(m)
    A[m:, :m] = D
    
    # Define semi-discrete approximation
    def rhs(u):
        return A@u
    
    # Define initial function
    def phi_0(x, y):
        return csc_matrix(np.exp(-(x-x_0)**2/o**2-(y-y_0)**2/o**2)).reshape((m, 1), order="C")

    # Set initial values (u = [phi, phi_t]^T)
    u = vstack((phi_0(X_vec, Y_vec), csc_matrix((m, 1)))).toarray()
    t = 0
    print(u)

    def SPL(u):
        return 20*np.log10(u/0.000002)


    # =====================
    # Plot with lab1 method
    # =====================
    # Initialize plot
    # fig, ax = plt.subplots()
    # zlow = -0.2
    # zhigh = 0.2
    # u_plot = u[:len(u)//2].reshape((m_y+1, m_x+1))
    # srf = ax.pcolor(X_vec, Y_vec, u_plot, shading = "nearest")#, vmin = zlow, vmax = zhigh)
    # # ax.set_xlim([x_l, x_r])
    # # ax.set_ylim([, ])
    # # fig.colorbar(srf, ax)
    # title = plt.title(f'Sound Pressure Level at Time T = {t:.3e}')
    # fig.tight_layout()
    # plt.draw()
    # plt.pause(0.5)

    # # Using RK4 to solve and store the values in phi_final_vec
    # for time_step in range(m_t):
    #     u, t = rk4.step(rhs, u, t, h_t)

    #     if show_time_steps and (time_step % 2 == 0 or time_step == m_t-1):
    #         u_plot = u[:len(u)//2].reshape((m_y+1, m_x+1))
    #         srf.remove()
    #         srf = ax.pcolor(X_vec, Y_vec, u_plot, shading = "nearest", vmin = zlow, vmax = zhigh)
    #         title.set_text(f'Sound Pressure Level at Time T = {t:.3e}')
    #         plt.draw()
    #         plt.pause(1e-6)
    

    # ========================================
    # Plot by saving in separate u_plot vector
    # ========================================
    # start_time = time.time()
    # u_plot = []

    # # Using RK4 to solve and store the values in phi_final_vec
    # for time_step in range(m_t):
    #     u, t = rk4.step(rhs, u, t, h_t)

    #     if show_time_steps:
    #         u_plot.append(u[:len(u)//2].reshape((m_y+1, m_x+1)))

    # if show_time_steps:
    #     fig, ax = plt.subplots()
    #     fig.colorbar(mappable = plt.scatter(np.zeros(10), np.linspace(0,1,10), s=30, c = np.linspace(0,1,10), 
    #                                             cmap='viridis'), ax = ax)
    #     for i, uplot in enumerate(u_plot):
    #         ax.clear()
    #         ax.contourf(X_vec, Y_vec, uplot, levels=30)
    #         ax.set_title(f'PDE at time T = {i*h_t}')
    #         ax.set_aspect('equal')
    #         plt.pause(1e-6)
    #     plt.show()
    # ===========================================
    # ===========================================
    
    # =============
    # Plot directly
    # =============
    # start_time = time.time()
    # if show_time_steps:
    #     fig, ax = plt.subplots()
    #     fig.colorbar(mappable = plt.scatter(np.zeros(10), np.linspace(0,1,10), s=30, c = np.linspace(0,1,10), 
    #                                             cmap='viridis'), ax = ax)

    # # Using RK4 to solve and store the values in phi_final_vec
    # for time_step in range(m_t):
    #     u, t = rk4.step(rhs, u, t, h_t)

    #     if show_time_steps:
    #         ax.clear()
    #         ax.contourf(X_vec, Y_vec, -u[len(u)//2:].reshape((m_y+1, m_x+1)), levels=30)
    #         ax.set_title(f'PDE at time T = {time_step*h_t}')
    #         ax.set_aspect('equal')
    #         plt.pause(1e-3)
    # plt.show()
    # ===========================================
    # ===========================================


    # ===============
    # Plot 3D surface
    # ===============
    # start_time = time.time()
    # if show_time_steps:
    #     ax = plt.axes(projection = '3d')

    # # Using RK4 to solve and store the values in phi_final_vec
    # for time_step in range(m_t):
    #     u, t = rk4.step(rhs, u, t, h_t)

    #     if show_time_steps and (time_step % 5 == 0 or time_step == m_t-1):
    #         ax.clear()
    #         ax.plot_surface(X_vec, Y_vec, -u[len(u)//2:].reshape((m_y+1, m_x+1)), rstride=1, cstride=1,
    #             cmap='viridis', edgecolor='none')
    #         ax.set_title(f'PDE at time T = {time_step*h_t:.3f}')
    #         ax.set_zlim([-15, 15])
    #         plt.pause(1e-6)
    # plt.show()
    # ===========================================
    # ===========================================

    # ===============
    # Plot 3D surface, alternative
    # ===============
    start_time = time.time()
    if show_time_steps:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection = '3d')

    # Using RK4 to solve and store the values in phi_final_vec
    for time_step in range(m_t):
        u, t = rk4.step(rhs, u, t, h_t)

        if show_time_steps and (time_step % 5 == 0 or time_step == m_t-1):
            ax.clear()
            ax.plot_surface(X_vec, Y_vec, -u[len(u)//2:].reshape((m_y+1, m_x+1)), rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
            ax.set_title(f'PDE at time T = {time_step*h_t:.3f}')
            ax.set_zlim([-15, 15])
            plt.pause(1e-6)
    plt.show()
    # ===========================================
    # ===========================================



    # ===========================
    # Plot 3D surface with mayavi
    # ===========================
    # start_time = time.time()

    # import sys
    # from vispy import app, scene, color
    # import matplotlib.colors

    # canvas = scene.SceneCanvas(keys='interactive', bgcolor='w')
    # view = canvas.central_widget.add_view()
    # view.camera = scene.TurntableCamera(up='z', fov=60)

    # # print(u[:len(u)//2].reshape((m_y+1, m_x+1)))
    # # print(-u[len(u)//2:].reshape((m_y+1, m_x+1)))

    # # U = -u[len(u)//2:].reshape((m_y+1, m_x+1))
    # # p1 = scene.visuals.SurfacePlot(X_vec, Y_vec, U, color=(0.3, 0.3, 1, 1))

    # # p1.transform = scene.transforms.MatrixTransform()
    # # p1.transform.scale([1/249., 1/249., 1/249.])
    # # p1.transform.translate([-0.5, -0.5, 0])

    # # view.add(p1)

    # # # Using RK4 to solve and store the values in phi_final_vec
    # for time_step in range(m_t):
    #     u, t = rk4.step(rhs, u, t, h_t)

    #     # if show_time_steps and (time_step % 5 == 0 or time_step == m_t-1):
    #         # p1 = scene.visuals.SurfacePlot(X_vec, Y_vec, -u[len(u)//2:].reshape((m_y+1, m_x+1)), color=(0.3, 0.3, 1, 1))
    #         # p1._update_data()
    #         # p1.set_data(z = U)

    

    # U = -u[len(u)//2:].reshape((m_y+1, m_x+1))
    # p1 = scene.visuals.SurfacePlot(X_vec, Y_vec, U, color=(1, 0, 0, 1))

    # p1.transform = scene.transforms.MatrixTransform()
    # p1.transform.scale([1, 1, 1/15.])
    # p1.transform.translate([0, 0, 0])

    # # p1.cmap = 'coolwarm'
    # # p1.clim = (U.min()/15, U.max()/15)

    # view.add(p1)

    # # cbar = scene.ColorBar(size=(100, 20), pos=(50, 10), cmap='coolwarm', label='Z')
    # # view.add(cbar)  

    # xax = scene.Axis(pos=[[-0.5, -0.5], [0.5, -0.5]], tick_direction=(0, -1),
    #                 font_size=16, axis_color='k', tick_color='k', text_color='k',
    #                 parent=view.scene)
    # xax.transform = scene.STTransform(translate=(0, 0, -0.2))

    # yax = scene.Axis(pos=[[-0.5, -0.5], [-0.5, 0.5]], tick_direction=(-1, 0),
    #                 font_size=16, axis_color='k', tick_color='k', text_color='k',
    #                 parent=view.scene)
    # yax.transform = scene.STTransform(translate=(0, 0, -0.2))

    # # Add a 3D axis to keep us oriented
    # # axis = scene.visuals.XYZAxis(parent=view.scene)

    # canvas.show()
    # if sys.flags.interactive == 0:
    #     app.run()
    # ===========================================
    # ===========================================


    time_ = time.time() - start_time
    print(time_)

    # # Reshape u and print final values
    # levels = 30
    # # levels = np.linspace(0, 140, 141)
    # u = -u[len(u)//2:].reshape((m_y+1, m_x+1))
    # plt.contourf(X_vec, Y_vec, u, levels = levels)
    # plt.title(f'Sound Pressure Level at Time T = {t:.3e}')
    # plt.xlabel('x')
    # plt.ylabel('y')
    # plt.gca().set_aspect('equal')
    # plt.colorbar()
    # plt.show()


if __name__ == '__main__':
    main()