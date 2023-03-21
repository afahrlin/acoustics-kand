"""4th order Runge-Kutta time-stepping.
This module solves ODEs of the form dv/dt = f(v).
"""

def step(f, v, t, dt):
    """Take one RK4 step. Return updated solution and time.
    f: Right-hand-side function: dv/dt = f(v)
    v: current solution
    t: current time
    dt: time step
    """

    # Compute rates k1-k4
    k1 = dt*f(v)
    k2 = dt*f(v + 0.5*k1)
    k3 = dt*f(v + 0.5*k2)
    k4 = dt*f(v + k3)

    # Update solution and time
    v = v + 1/6*(k1 + 2*k2 + 2*k3 + k4)
    t = t + dt

    return v, t