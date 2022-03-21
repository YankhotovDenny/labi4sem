from fenics import *
import os
import numpy as np

T = 200.0            # final time
num_steps = 50     # number of time steps
dt = T / num_steps # time step size
alpha1 = 1
alpha2 = 1          # parameter alpha
beta = 1         # parameter beta

# Create mesh and define function space
mesh = Mesh("iko3.xml")
V = FunctionSpace(mesh, 'P', 1)

# Define boundary condition
u_D = Expression('abs(x[0])*x[0] + alpha1*abs(x[1])*x[1] + alpha2*abs(x[2])*x[2] + beta*t',
                 degree=3, alpha1=alpha1, alpha2=alpha2, beta=beta, t=0)

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, u_D, boundary)

# Define initial value
u_n = interpolate(u_D, V)
#u_n = project(u_D, V)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(beta - 2 - 2*(alpha1 + alpha2))

F = u*v*dx + dt*dot(grad(u), grad(v))*dx - (u_n + dt*f)*v*dx
a, L = lhs(F), rhs(F)

# Time-stepping
u = Function(V)
t = 0

res_file = File('heat/solution.pvd')

for n in range(num_steps):

    # Update current time
    t += dt
    u_D.t = t

    # Compute solution
    solve(a == L, u, bc)

    # Plot solution
    plot(u)

    # Compute error at vertices
    u_e = interpolate(u_D, V)
    error = np.abs(u_e.vector() - u.vector()).max()
    print('t = %.2f: error = %.3g' % (t, error))
    
    res_file << u

    # Update previous solution
    u_n.assign(u)
