import numpy as np
import matplotlib.pyplot as plt
import multiprocessing

# Domain definition
a = 0
b = 1

# Domain discretization
npts = 200
x = np.linspace(a,b,npts)

# Initial condition
u = np.zeros_like(x)
u_mask = (x > 0.4) & (x < 0.7)
u[u_mask] = 1.0

# Advection info
a = -1.0

# Solver
dt = 0.01
dh = abs(a-b)/npts
t = 0
tf = 10

# CFL condition 0 <= a*dt/dh <= 1 
CFL = 0.1
dt = CFL*dh/abs(a)
print("a*dt/dh = ", abs(a)*dt/dh)
print("dt = ", dt)

# Create the plot
plt.ion()  # Turn on interactive mode for live updating
fig, ax = plt.subplots()
line, = ax.plot(x, u)

u_next = np.zeros_like(u)

def f(u):
    return np.multiply(u,u)

while t<tf:
    print(t)
    
    if a>0:
        for i in range(1,u.shape[0]):
            u_next[i] = u[i] - a*(dt/dh)*(f(u[i]) - f(u[i-1]))
        u_next[0] = u[-1]

    if a<0:
        for i in range(0,u.shape[0]-1):
            u_next[i] = u[i] - a*(dt/dh)*(f(u[i+1]) - f(u[i]))
        u_next[-1] = u[0]

    t+=dt
    u = u_next

    line.set_ydata(u_next)  # Update the y-data of the line
    plt.draw()  # Redraw the plot
    plt.pause(0.01)  # Add a small delay to control the update rate

plt.ioff()  # Turn off interactive mode when done
plt.show()  # Display the final plot