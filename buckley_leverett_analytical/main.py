from relative_permeability import Relative_perm
from fractional_flow import Fractional_flow
import numpy as np
import matplotlib.pyplot as plt

# Wettability scenarios

# Mixed-wet
kr_infos = {'Swc' : 0.1,
            'Sor' : 0.15,
            'krw_max': 0.5,
            'kro_max': 1.0,
            'a': 8,
            'b': 2.5}

# # Strongly water-wet
# kr_infos = {'Swc' : 0.1,
#             'Sor' : 0.4,
#             'krw_max': 0.1,
#             'kro_max': 1.0,
#             'a': 2,
#             'b': 1}

# Define relative permeability functions
rel_perm = Relative_perm(kr_infos)

# Define Sw range
Sw = np.linspace(kr_infos['Swc'],1-kr_infos['Sor'],200)

# Builds and show relative permeability vector values
rel_perm.set_krw(Sw)
rel_perm.set_kro(Sw)
rel_perm.plot_perm(Sw)

# Definig fractional flow
fw_infos = {'M' : 1.0}  # M = mu_displaced / mu_injected
frac_flow = Fractional_flow(fw_infos)
frac_flow.set_frac_flow(Sw, rel_perm.get_krw(), rel_perm.get_kro())
frac_flow.plot_fw(Sw)


