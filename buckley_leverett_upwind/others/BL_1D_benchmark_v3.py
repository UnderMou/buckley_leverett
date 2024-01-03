import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def f(u,c):
    return u*u/(u*u + c*(1 - u)*(1 - u))

def dfdu(u,c):
    return (2*c*(1-u)*u) / (np.power( c*(u-1)*(u-1) + u*u ,2))

data = pd.read_csv("bl_upwind_t0_5.csv", delimiter=',', index_col = False)

# Domain definition
xi = 0.0
xf = 1.0

# Domain discretization
nel = 200
x = np.linspace(xi,xf,nel+1)
dx = abs(xi-xf)/(nel)

# Time discretization
ti = 0.0
tf = 0.5
dt = 0.002
dfdu_max = 2.0
print("CFLmax = ", dfdu_max*dt/dx)

# Solution initialization
u = np.zeros_like(x)    # Initial condition
u[0] = 1.0              # Boundary condition at x = 0 | u = 1.0

# Create the plot
plt.ion()  # Turn on interactive mode for live updating
fig, ax = plt.subplots()
ax.set_ylim([-0.1, 1.1])
line, = ax.plot(x, u)

# Fractional flux function 
c = 1.0

# Solver
u_next = np.copy(u)
t = ti
while t<=tf:
    print(t)
    
    u_next[0] = np.copy(u[0])   # DIRICHLET B.C.

    ######################################################

    # for i in range(1,u.shape[0]-1):
        
    #     if u[i] != u[i+1]: v_next = (dt/dx)*(f(u[i+1],c) - f(u[i],c))/(u[i+1]-u[i])
    #     else: v_next = (dt/dx)*dfdu(u[i],c)

    #     if u[i] != u[i-1]: v_prev = (dt/dx)*(f(u[i],c) - f(u[i-1],c))/(u[i]-u[i-1])
    #     else: v_prev = (dt/dx)*dfdu(u[i],c)

    #     if v_next > 0: H_next = f(u[i],c)
    #     else: H_next = f(u[i+1],c)

    #     if v_prev > 0: H_prev = f(u[i-1],c)
    #     else: H_prev = f(u[i],c)

    #     u_next[i] = u[i] - (dt/dx)*(H_next - H_prev)

    ######################################################

    a = 1.0
    for i in range(1,u.shape[0]):
        u_next[i] = u[i] - a*(dt/dx)*(f(u[i],c) - f(u[i-1],c))

    ######################################################

    u_next[-1] = u_next[-2]     # NEUMANN B.C.  

    line.set_ydata(u)   # Update the y-data of the line
    plt.draw()          # Redraw the plot
    plt.pause(0.01)     # Add a small delay to control the update rate

    t+=dt
    u = u_next.copy()
    # u = u_next

fig, ax2 = plt.subplots()
ax2.plot(x,u,c='b',label='This code')
ax2.plot(data['x'],data['u'],c='r',label='Reference')
plt.grid(True)
plt.legend()
plt.show()

plt.ioff()  # Turn off interactive mode when done
plt.show()  # Display the final plot