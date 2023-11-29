import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from relative_permeability import Relative_perm
from fractional_flow import Fractional_flow

# def f(u,c):
#     return u*u/(u*u + c*(1 - u)*(1 - u))

# def dfdu(u,c):
#     return (2*c*(1-u)*u) / (np.power( c*(u-1)*(u-1) + u*u ,2))

# data = pd.read_csv("bl_upwind_t0_5.csv", delimiter=',', index_col = False)



# Wettability scenarios

# # Strongly water-wet
# kr_infos = {'Swc' : 0.1,
#             'Sor' : 0.4,
#             'krw_max': 0.1,
#             'kro_max': 1.0,
#             'a': 2,
#             'b': 1,
#             'M' : 200.0,
#             'wettability' : 'Strongly-water-wet'}

# Mixed-wet
kr_infos = {'Swc' : 0.1,
            'Sor' : 0.15,
            'krw_max': 0.5,
            'kro_max': 1.0,
            'a': 8,
            'b': 2.5,
            'M' : 200.0,
            'wettability' : 'Mixed-wet'}

# Define Sw range
Sw = np.linspace(kr_infos['Swc'],1-kr_infos['Sor'],2000)

# Define relative permeability functions
rel_perm = Relative_perm(Sw, kr_infos)
# Builds and show relative permeability vector values
rel_perm.set_krw()
rel_perm.set_kro()
rel_perm.plot_perm(save_pdf = False)

# Definig fractional flow
frac_flow = Fractional_flow(Sw, rel_perm, kr_infos)
frac_flow.set_fw()
frac_flow.plot_fw(save_pdf = False)


# Domain definition
xi = 0.0
xf = 1.0

# Domain discretization
nel = 300
x = np.linspace(xi,xf,nel+1)
dx = abs(xi-xf)/(nel)

# Time discretization
ti = 0.0
tf = 8000
dt = 1.0
dfdu_max = 2.0
print("CFLmax = ", dfdu_max*dt/dx)

# Solution initialization
u = np.ones_like(x)*kr_infos['Swc']                 # Initial condition
u[0] = 1-kr_infos['Sor']              # Boundary condition at x = 0 | u = 1.0

# Create the plot
plt.ion()  # Turn on interactive mode for live updating
fig, ax = plt.subplots()
ax.set_ylim([-0.1, 1.1])
ax.grid(True)
ax.set_xlabel('x')
ax.set_ylabel('Sw')
line, = ax.plot(x, u)

text_element = ax.text(0.825, 1.05, '', fontsize=10, ha='left', va='center', color='k')

dimensional_reservoir = {'L'  : 1.0,
                         'phi': 0.35,
                         'qt' : 9.5e-5,
                         'ti' : 0.01,
                         'tf' : 7000,
                         'Nt' : 200,
                         'Nx' : 1000,
                         'A'  : 1.0} 

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

    for i in range(1,u.shape[0]):
        u_next[i] = u[i] - (dimensional_reservoir['qt']/dimensional_reservoir['phi'])*(dt/dx)*(frac_flow.eval_fw(u[i]) - frac_flow.eval_fw(u[i-1]))

    ######################################################  
    
    new_text = f'time: {t:.2f}'
    text_element.set_text(new_text)

    line.set_ydata(u)   # Update the y-data of the line
    plt.draw()          # Redraw the plot
    plt.pause(0.01)     # Add a small delay to control the update rate

    t+=dt
    u = u_next.copy()
    # u = u_next

fig, ax2 = plt.subplots()
ax2.plot(x,u,c='b',label='This code')
# ax2.plot(data['x'],data['u'],c='r',label='Reference')
plt.grid(True)
plt.legend()
plt.show()

plt.ioff()  # Turn off interactive mode when done
plt.show()  # Display the final plot