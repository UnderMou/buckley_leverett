import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv("bl_upwind_t0_5.csv", delimiter=',', index_col = False)

# Domain definition
xi = 0
xf = 1

# Domain discretization
npts = 50
x = np.linspace(xi,xf,npts)
print(x)

# Initial condition
u = np.zeros_like(x)
u[0] = 1.0


# Solver
dh = abs(xi-xf)/(npts-1)
print("dh = ", dh)

# CFL condition 0 <= a*dt/dh <= 1 
t = 0.0
tf = 10

dt = 0.008
# CFL = 0.1
# dt = CFL*dh/abs(2.0)
# print("a*dt/dh = ", abs(2.0)*dt/dh)
# print("dt = ", dt)

# Create the plot
plt.ion()  # Turn on interactive mode for live updating
fig, ax = plt.subplots()
ax.set_ylim([-0.1, 1.1])
line, = ax.plot(x, u)

u_next = np.copy(u)

# buckley leverett fractional flux
c = 1.0
def f(u):
    return u*u/(u*u + c*(1 - u)*(1 - u))

def dfdu(u):
    #return (2*u*(1-u))/((2*u*u - 2*u +1)*(2*u*u - 2*u +1))
    return (2*c*(1-u)*u) / (np.power( c*(u-1)*(u-1) + u*u ,2))

while t<tf:
    print(t)
    
    u_next[0] = 1.0

    for i in range(1,u.shape[0]-1):
        
        if u[i] != u[i+1]: v_next = (dt/dh)*(f(u[i+1]) - f(u[i]))/(u[i+1]-u[i])
        else: v_next = (dt/dh)*dfdu(u[i])

        if u[i] != u[i-1]: v_prev = (dt/dh)*(f(u[i]) - f(u[i-1]))/(u[i]-u[i-1])
        else: v_prev = (dt/dh)*dfdu(u[i])

        if v_next > 0: H_next = f(u[i])
        else: H_next = f(u[i+1])

        if v_prev > 0: H_prev = f(u[i-1])
        else: H_prev = f(u[i])

        u_next[i] = u[i] - (dt/dh)*(H_next - H_prev)
    
    u_next[-1] = u_next[-2]

    line.set_ydata(u)  # Update the y-data of the line
    plt.draw()  # Redraw the plot
    plt.pause(0.01)  # Add a small delay to control the update rate

    if (t > 0.5) :
        fig, ax2 = plt.subplots()
        ax2.plot(x,u,c='b')
        ax2.plot(data['x'],data['u'],c='r')
        plt.show()
        break

    t+=dt
    u = u_next

plt.ioff()  # Turn off interactive mode when done
plt.show()  # Display the final plot