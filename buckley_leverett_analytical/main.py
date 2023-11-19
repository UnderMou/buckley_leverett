from relative_permeability import Relative_perm
from fractional_flow import Fractional_flow
from bl_solution import BuckleyLeverettSolution

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
#             'M' : 200.0,
#             'wettability' : 'Mixed-wet'}

# Strongly water-wet
kr_infos = {'Swc' : 0.1,
            'Sor' : 0.4,
            'krw_max': 0.1,
            'kro_max': 1.0,
            'a': 2,
            'b': 1,
            'M' : 200.0,
            'wettability' : 'Mixed-wet'}

# # Weakly water-wet
# kr_infos = {'Swc' : 0.1,
#             'Sor' : 0.3,
#             'krw_max': 0.2,
#             'kro_max': 1.0,
#             'a': 2,
#             'b': 1.5,
#             'M' : 200.0,
#             'wettability' : 'Weakly water-wet'}

# # Oil-wet
# kr_infos = {'Swc' : 0.1,
#             'Sor' : 0.05,
#             'krw_max': 0.95,
#             'kro_max': 1.0,
#             'a': 1.5,
#             'b': 4,
#             'M' : 200.0,
#             'wettability' : 'Oil-wet'}

###########################################################################################

# Define Sw range
Sw = np.linspace(kr_infos['Swc'],1-kr_infos['Sor'],2000)

# Define relative permeability functions
rel_perm = Relative_perm(Sw, kr_infos)
# Builds and show relative permeability vector values
rel_perm.set_krw()
rel_perm.set_kro()
rel_perm.plot_perm()

# Definig fractional flow
frac_flow = Fractional_flow(Sw, rel_perm, kr_infos)
frac_flow.set_fw()
frac_flow.plot_fw()

# Construction of analytical solution
bl_solution = BuckleyLeverettSolution(Sw, frac_flow, kr_infos)
bl_solution.construct_solution()
bl_solution.show_solution()












# # Recovery Calculations

# Sw1 = np.linspace(Sws, 1-kr_infos['Sor'], 100)
# Sw1 = Sw1[1:len(Sw1)-1]

# NpD = np.zeros_like(Sw1)
# tD1 = np.zeros_like(Sw1)

# for i in range(len(Sw1)):
#     Sw1_idx = np.argmin(np.abs(Sw - Sw1[i]))
#     # print(Sw[Sw1_idx], Sw1[i])
#     fw1 = frac_flow.get_fw()[Sw1_idx]
#     # print(Sw1[i], fw1)
#     tD1[i] = 1 / dfw_dSw[Sw1_idx]
#     # print(tD1[i])

#     Sw_bar = Sw1[i] + tD1[i]*(1 - fw1)

#     NpD[i] = Sw_bar - kr_infos['Swc']

# tDmax = 3.0

# id_tDmax = np.argmin(np.abs(tDmax - tD1))

# tD = np.concatenate((np.array([0.0, 1/vsD]), tD1[:id_tDmax])) 
# NpD = np.concatenate((np.array([0.0, 1/vsD]), NpD[:id_tDmax])) 

# plt.plot(tD,NpD,c='b',label='Analítico ' + r'$\mu_{disp.}/\mu_{injec.} = $' + str(kr_infos['M']))
# plt.ylim([-0.1, 1.1])
# plt.grid(True)
# plt.title('Wettability condition: ' + kr_infos['wettability'])
# plt.ylabel(r"$N_{p_D}$")
# plt.xlabel(r"$t_D$")
# plt.legend()
# plt.show()