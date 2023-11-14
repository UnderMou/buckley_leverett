from relative_permeability import Relative_perm
from fractional_flow import Fractional_flow
import numpy as np
import matplotlib.pyplot as plt

# Wettability scenarios

# # Mixed-wet
# kr_infos = {'Swc' : 0.1,
#             'Sor' : 0.15,
#             'krw_max': 0.5,
#             'kro_max': 1.0,
#             'a': 8,
#             'b': 2.5}

# # Strongly water-wet
# kr_infos = {'Swc' : 0.1,
#             'Sor' : 0.4,
#             'krw_max': 0.1,
#             'kro_max': 1.0,
#             'a': 2,
#             'b': 1}

# # Weakly water-wet
# kr_infos = {'Swc' : 0.1,
#             'Sor' : 0.3,
#             'krw_max': 0.2,
#             'kro_max': 1.0,
#             'a': 2,
#             'b': 1.5}

# Weakly water-wet
kr_infos = {'Swc' : 0.1,
            'Sor' : 0.05,
            'krw_max': 0.95,
            'kro_max': 1.0,
            'a': 1.5,
            'b': 4}

# Define relative permeability functions
rel_perm = Relative_perm(kr_infos)

# Define Sw range
Sw = np.linspace(kr_infos['Swc'],1-kr_infos['Sor'],2000)

# Builds and show relative permeability vector values
rel_perm.set_krw(Sw)
rel_perm.set_kro(Sw)
# rel_perm.plot_perm(Sw)

# Definig fractional flow
fw_infos = {'M' : 1.0}  # M = mu_displaced / mu_injected
frac_flow = Fractional_flow(fw_infos)
frac_flow.set_fw(Sw, rel_perm.get_krw(), rel_perm.get_kro())
# frac_flow.plot_fw(Sw)

# Compute the derivative using NumPy's gradient function
dfw_dSw = np.gradient(frac_flow.get_fw(), Sw)
# plt.plot(Sw,dfw_dSw)
# plt.show()



# Constructing the solution

# Compute Sws and shock solution
# TODO: Pensar numa abordagem melhor. Está muito gambiarra!
diff = np.abs(np.divide(frac_flow.get_fw(), Sw - kr_infos['Swc']+1e-6) - dfw_dSw) + 1-Sw
# plt.plot(Sw,diff)
# plt.show()
Sws_idx = np.argmin(diff)
Sws = Sw[Sws_idx]
# print(Sws)

# Compute vsD
vsD = frac_flow.get_fw()[Sws_idx] / (Sws - kr_infos['Swc'])
# print(vsD)



# Compute solution constant state at vD <= vD_min
vD_min = dfw_dSw[-1]

vD1 = np.linspace(0,vD_min,100)
Sw1 = 1-kr_infos['Sor'] * np.ones_like(vD1)
# print(Sw1)
# print(vD1)



# Compute rarefaction solution
Sw2 = Sw[Sws_idx:][::-1]
Sw2[-1] = Sws
# print(Sw2)
# Sw2 = np.concatenate((Sw2,np.array([Sws])))
# print(len(Sw_rare))
vD2 = dfw_dSw[Sws_idx:][::-1]
vD2[-1] = vsD
# print(vD2)
# vD2 = np.concatenate((vD2,np.array([vsD])))
# print(len(vD_rare))
# print(vD2,Sw2)



# Compute solution constant state at vD > vsD
vD3 = np.linspace(vsD,6.0,500)
Sw3 = kr_infos['Swc'] * np.ones_like(vD3)



# Whole solution grouping
Sol_vD = np.concatenate((vD1 , vD2, vD3))
Sol_Sw = np.concatenate((Sw1 , Sw2, Sw3))

plt.plot(Sol_vD,Sol_Sw)
plt.ylim([-0.1,1.1])
plt.xlim([-0.1,6.0])
plt.grid(True)
plt.show()