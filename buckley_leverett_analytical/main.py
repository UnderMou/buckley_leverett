from relative_permeability import Relative_perm
from fractional_flow import Fractional_flow
from bl_solution import BuckleyLeverettSolution
from recovery_calc import Recovery_calc

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

# Strongly water-wet
kr_infos = {'Swc' : 0.1,
            'Sor' : 0.4,
            'krw_max': 0.1,
            'kro_max': 1.0,
            'a': 2,
            'b': 1,
            'M' : 20.0,
            'wettability' : 'Strongly-water-wet'}

# # Weakly water-wet
# kr_infos = {'Swc' : 0.1,
#             'Sor' : 0.3,
#             'krw_max': 0.2,
#             'kro_max': 1.0,
#             'a': 2,
#             'b': 1.5,
#             'M' : 20.0,
#             'wettability' : 'Weakly-water-wet'}

# # Oil-wet
# kr_infos = {'Swc' : 0.1,
#             'Sor' : 0.05,
#             'krw_max': 0.95,
#             'kro_max': 1.0,
#             'a': 1.5,
#             'b': 4,
#             'M' : 20.0,
#             'wettability' : 'Oil-wet'}

###########################################################################################

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

# Construction of recovery curve
rec_curve = Recovery_calc(kr_infos, frac_flow, bl_solution)
rec_curve.do_recovery()
# rec_curve.show_curve(save_pdf = True)

# Recovery of dimensional solution
dimensional_reservoir = {'L'  : 1.0,
                         'phi': 0.35,
                         'qt' : 1e-5,
                         'ti' : 0.01,
                         'tf' : 40000,
                         'Nt' : 500,
                         'Nx' : 1000}

bl_solution.do_dimensional_Sw_x(dimensional_reservoir)
bl_solution.show_dimensional_Sw_x()

rec_curve.do_dimensional_NpD_t(dimensional_reservoir)
rec_curve.show_dimensional_NpD_t()