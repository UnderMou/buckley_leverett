import numpy as np
import matplotlib.pyplot as plt

# Domain definition
xi = 0
xf = 1

# Domain discretization
npts = 50
x = np.linspace(xi,xf,npts)
dh = abs(xi-xf)/npts
u = np.zeros_like(x)
u[0] = 1.0

# Create the plot
plt.ion()  # Turn on interactive mode for live updating
fig, ax = plt.subplots()
ax.set_ylim([-0.1, 1.1])
line, = ax.plot(x, u)

# Solution
t = 0.1
tf = 10
dt = 0.1

while t<tf:
    print(t)
    
    for i in range(1,u.shape[0]):
        if x[i]/t < 0.5*(1+np.sqrt(2)):
            A = -2*x[i]/t
            B = np.sqrt(4*x[i]/t+1)
            C = x[i]/t
            u[i] = 0.5*(np.sqrt( (A+B-1) / (C) + 1) +1) 

    line.set_ydata(u)  # Update the y-data of the line
    plt.draw()  # Redraw the plot
    plt.pause(0.01)  # Add a small delay to control the update rate

    if t >= 0.5:
        break

    t+=dt

plt.ioff()  # Turn off interactive mode when done
plt.show()  # Display the final plot








u = np.zeros_like(x)
# u_mask = (x > 0.4) & (x < 0.7)
# u[u_mask] = 1.0

# Advection info
a = 1.0

# Solver
dt = 0.01
dh = abs(xi-xf)/npts
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
ax.set_ylim([-0.1, 1.1])
line, = ax.plot(x, u)

u_next = np.zeros_like(u)

# buckley leverett fractional flux
c = 1.0

def f(u):
    return u*u/(u*u + c*np.power(1 - u, 2))

while t<tf:
    print(t)
    
    if a>0:
        u_next[0] = 1.0
        for i in range(1,u.shape[0]):
            u_next[i] = u[i] - a*(dt/dh)*(f(u[i]) - f(u[i-1]))
        

    if a<0:
        u_next[u.shape[0]-1] = 1.0
        for i in range(0,u.shape[0]-1):
            u_next[i] = u[i] - a*(dt/dh)*(f(u[i+1]) - f(u[i]))

    t+=dt
    u = u_next

    line.set_ydata(u_next)  # Update the y-data of the line
    plt.draw()  # Redraw the plot
    plt.pause(0.01)  # Add a small delay to control the update rate

plt.ioff()  # Turn off interactive mode when done
plt.show()  # Display the final plot