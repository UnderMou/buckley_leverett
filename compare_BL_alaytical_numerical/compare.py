from analytical_model.relative_permeability import Relative_perm
from analytical_model.fractional_flow import Fractional_flow
from analytical_model.bl_solution import BuckleyLeverettSolution
from analytical_model.recovery_calc import Recovery_calc

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

###########################################################################################

# Wettability scenarios

# # Mixed-wet
# kr_infos = {'Swc' : 0.1,
#             'Sor' : 0.15,
#             'krw_max': 0.5,
#             'kro_max': 1.0,
#             'a': 8,
#             'b': 2.5,
#             'M' : 20.0,
#             'wettability' : 'Mixed-wet'}

# # Strongly water-wet
# kr_infos = {'Swc' : 0.1,
#             'Sor' : 0.4,
#             'krw_max': 0.1,
#             'kro_max': 1.0,
#             'a': 2,
#             'b': 1,
#             'M' : 200.0,
#             'wettability' : 'Strongly-water-wet'}

# Weakly water-wet
kr_infos = {'Swc' : 0.1,
            'Sor' : 0.3,
            'krw_max': 0.2,
            'kro_max': 1.0,
            'a': 2,
            'b': 1.5,
            'M' : 200.0,
            'wettability' : 'Weakly-water-wet'}

# # Oil-wet
# kr_infos = {'Swc' : 0.1,
#             'Sor' : 0.05,
#             'krw_max': 0.95,
#             'kro_max': 1.0,
#             'a': 1.5,
#             'b': 4,
#             'M' : 0.005,
#             'wettability' : 'Oil-wet'}

###########################################################################################

# Recovery of dimensional solutions
dimensional_reservoir = {'L'  : 1.0,
                         'phi': 0.35,
                         'qt' : 7.5e-5,
                         'ti' : 0.01,
                         'tf' : 500,
                         'Nt' : 200,
                         'Nx' : 1000,
                         'A'  : 1.0} 

# ANALYTICAL MODEL

print("doing analytical ...")
# Define Sw range
Sw = np.linspace(kr_infos['Swc'],1-kr_infos['Sor'],2000)
# Define relative permeability functions
rel_perm = Relative_perm(Sw, kr_infos)
# Builds and show relative permeability vector values
rel_perm.set_krw()
rel_perm.set_kro()
# rel_perm.plot_perm(save_pdf = True)
# Definig fractional flow
frac_flow = Fractional_flow(Sw, rel_perm, kr_infos)
frac_flow.set_fw()
# frac_flow.plot_fw(save_pdf = True)
# Construction of analytical solution
bl_solution = BuckleyLeverettSolution(Sw, frac_flow, kr_infos)
bl_solution.construct_solution()
# bl_solution.show_solution(save_pdf = True)
# # Construction of recovery curve
# rec_curve = Recovery_calc(kr_infos, frac_flow, bl_solution)
# rec_curve.do_recovery()
# # rec_curve.show_curve(save_pdf = False)
# Dimensional solution
bl_solution.do_dimensional_Sw_x(dimensional_reservoir)



# NUMERICAL MODEL

print("doing numerical ...")
# Domain definition
xi = 0.0
xf = dimensional_reservoir['L']
# Domain discretization
nel = 300
x = np.linspace(xi,xf,nel+1)
dx = abs(xi-xf)/(nel)
# Time discretization
ti = dimensional_reservoir['ti']
tf = dimensional_reservoir['tf']
dt = 1.0
# Initial condition and dirichlet boundary condition
u = np.ones_like(x)*kr_infos['Swc']     # Initial condition
u[0] = 1-kr_infos['Sor']                # Boundary condition at x = 0 | u = 1 - Sor
# Solver
u_next = np.copy(u)
t = ti
while t<=tf:
    
    u_next[0] = np.copy(u[0])   # DIRICHLET B.C.

    for i in range(1,u.shape[0]):
        u_next[i] = u[i] - (dimensional_reservoir['qt']/dimensional_reservoir['phi'])*(dt/dx)*(frac_flow.eval_fw(u[i]) - frac_flow.eval_fw(u[i-1]))

    t+=dt
    u = u_next.copy()


# SHOW RESULTS

plt.scatter(x,u,c='r',s=5,label='numerical')
plt.plot(bl_solution.get_x(), bl_solution.get_grid()[-1][:], c='b', label='analytical')
plt.legend()
plt.xlabel('x')
plt.ylabel('Sw')
plt.grid(True)
plt.show()