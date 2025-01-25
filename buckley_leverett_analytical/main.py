from relative_permeability import Relative_perm
from fractional_flow import Fractional_flow, Fractional_flow_grav
from bl_solution import BuckleyLeverettSolution
from recovery_calc import Recovery_calc

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


###########################################################################################

# Wettability scenarios

# Mixed-wet
kr_infos = {'Swc' : 0,
            'Sor' : 0,
            'krw_max': 1,
            'kro_max': 1,
            'a': 2,
            'b': 2,
            'M' : 5,
            'wettability' : 'Durlofsky test case',
            'muw': 2.94e-14,
            'rho_w': 1.0,
            'rho_o': 0.7,
            'qt': 1.0,
            'K': 1e-13,
            'g': -9.81}

# # Strongly water-wet
# kr_infos = {'Swc' : 0.1,
#             'Sor' : 0.4,
#             'krw_max': 0.1,
#             'kro_max': 1.0,
#             'a': 2,
#             'b': 1,
#             'M' : 20.0,
#             'wettability' : 'Strongly-water-wet'}

# # Weakly water-wet
# kr_infos = {'Swc' : 0.1,
#             'Sor' : 0.3,
#             'krw_max': 0.2,
#             'kro_max': 1.0,
#             'a': 2,
#             'b': 1.5,
#             'M' : 200.0,
#             'wettability' : 'Weakly-water-wet'}

# # Oil-wet
# kr_infos = {'Swc' : 0.1,
#             'Sor' : 0.05,
#             'krw_max': 0.95,
#             'kro_max': 1.0,
#             'a': 1.5,
#             'b': 4,
#             'M' : 1.0,
#             'wettability' : 'Oil-wet'}

###########################################################################################

# Define Sw range
Sw = np.linspace(kr_infos['Swc'],1-kr_infos['Sor'],2000)

# Define relative permeability functions
rel_perm = Relative_perm(Sw, kr_infos)

# Builds and show relative permeability vector values
rel_perm.set_krw()
rel_perm.set_kro()
#rel_perm.plot_perm(save_pdf = False)

# Definig fractional flow
# frac_flow = Fractional_flow(Sw, rel_perm, kr_infos)
frac_flow = Fractional_flow_grav(Sw, rel_perm, kr_infos)
frac_flow.set_fw()
frac_flow.plot_fw(save_pdf = True)

# Construction of analytical solution
bl_solution = BuckleyLeverettSolution(Sw, frac_flow, kr_infos)
bl_solution.construct_solution()
bl_solution.show_solution(save_pdf = True)

# # Construction of recovery curve
# rec_curve = Recovery_calc(kr_infos, frac_flow, bl_solution)
# rec_curve.do_recovery()
# rec_curve.show_curve(save_pdf = False)

# tfind = 6250/1e5
tfind = 12500/1e5
# Recovery of dimensional solutions
dimensional_reservoir = {'L'  : 1.0,
                         'phi': 0.25,
                         'qt' : 1,
                         'ti' : 1e-5,
                         'tf' : tfind,
                         'Nt' : 800,
                         'Nx' : 1000,
                         'A'  : 1.0} 

bl_solution.do_dimensional_Sw_x(dimensional_reservoir)
bl_solution.show_dimensional_Sw_x()

# print(bl_solution.grid.shape)
# np.savetxt("Sw_analitycal_Grav_t025.csv", bl_solution.grid[-1,:], delimiter=",", fmt="%.4f")
# np.savetxt("analitycal_x.csv", bl_solution.x, delimiter=",", fmt="%.4f")

# rec_curve.do_dimensional_NpD_t(dimensional_reservoir)
# rec_curve.show_dimensional_NpD_t(dimensional_reservoir)

# rec_curve.define_integrand(bl_solution, frac_flow)
# rec_curve.production_integration(bl_solution, dimensional_reservoir)
# rec_curve.show_production_prod_t(dimensional_reservoir)

