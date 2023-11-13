from relative_permeability import Relative_perm
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
Sw = np.linspace(kr_infos['Swc'],1-kr_infos['Sor'],100)

# Builds and show relative permeability vector values
krw = rel_perm.set_krw(Sw)
kro = rel_perm.set_kro(Sw)
rel_perm.plot_perm(Sw)



