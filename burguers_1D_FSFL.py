from re import S
import numpy as np
import matplotlib.pyplot as plt
import multiprocessing

# Domain definition
x0 = 0
xf = 1

# Domain discretizations
npts = 100
x = np.linspace(x0,xf,npts)

# Initial condition
u = np.zeros_like(x)
u_mask = (x > 0.4) & (x < 0.7)
u[u_mask] = 1.0

# Advection info
a = -1.0

# Solver
dt = 0.01
dh = abs(xf-x0)/npts
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

def ufront(u,i):
    if i>=u.shape[0]: veloc = u[u.shape[0]-1]
    elif i<=0: veloc = u[0]
    else: veloc = u[i]

    return veloc

def FSFL(uU,uD,uR):
    beta = 2

    uUch = (uU-uR)/(uD-uR)

    if (uUch>=0) & (uUch<=1): ufc = uR + (uD-uR)*((-2*beta+4)*pow(uUch,4) + (4*beta-8)*pow(uUch,3) + 0.5*(-5*beta+8)*pow(uUch,2) + 0.5*(beta+2)*uUch)
    else : ufc = uU

    return ufc

def CONV(a,u,i):
    uf = (a + a)/2
    ug = (a + a)/2

    # computational face f
    if uf>0:
        uD = ufront(u,i+1)
        uR = ufront(u,i-1)
        uU = ufront(u,i)
        ufc = FSFL(uU,uD,uR)
    elif uf<0:
        uD = ufront(u,i)
        uR = ufront(u,i+2)
        uU = ufront(u,i+1)
        ufc = FSFL(uU,uD,uR)

    # computational face g
    if ug>0:
        uD = ufront(u,i)
        uR = ufront(u,i-2)
        uU = ufront(u,i-1)
        ugc = FSFL(uU,uD,uR)
    elif ug<0:
        uD = ufront(u,i-1)
        uR = ufront(u,i+1)
        uU = ufront(u,i)
        ugc = FSFL(uU,uD,uR)
    
    CONV = (uf*ufc - ug*ugc)
    return CONV

def f(u):
    return np.multiply(u,u)

u_next = np.zeros_like(u)
while t<tf:
    print(t)
    
    if a>0:
        for i in range(1,u.shape[0]):
            u_next[i] = u[i] - a*(dt/dh)*CONV(a,f(u),i)
        u_next[0] = u[-1]

    if a<0:
        for i in range(0,u.shape[0]-1):
            u_next[i] = u[i] + a*(dt/dh)*CONV(a,f(u),i)
        u_next[-1] = u[0]

    t+=dt
    u = u_next

    line.set_ydata(u_next)  # Update the y-data of the line
    plt.draw()  # Redraw the plot
    plt.pause(0.01)  # Add a small delay to control the update rate

plt.ioff()  # Turn off interactive mode when done
plt.show()  # Display the final plot