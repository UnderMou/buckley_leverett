from relative_permeability import Relative_perm
from fractional_flow import Fractional_flow
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Wettability scenarios

# Mixed-wet
kr_infos = {'Swc' : 0.1,
            'Sor' : 0.15,
            'krw_max': 0.5,
            'kro_max': 1.0,
            'a': 8,
            'b': 2.5,
            'wettability' : 'Mixed-wet'}

# # Strongly water-wet
# kr_infos = {'Swc' : 0.1,
#             'Sor' : 0.4,
#             'krw_max': 0.1,
#             'kro_max': 1.0,
#             'a': 2,
#             'b': 1,
#             'wettability' : 'Mixed-wet'}

# # Weakly water-wet
# kr_infos = {'Swc' : 0.1,
#             'Sor' : 0.3,
#             'krw_max': 0.2,
#             'kro_max': 1.0,
#             'a': 2,
#             'b': 1.5,
#             'wettability' : 'Weakly water-wet'}

# # Oil-wet
# kr_infos = {'Swc' : 0.1,
#             'Sor' : 0.05,
#             'krw_max': 0.95,
#             'kro_max': 1.0,
#             'a': 1.5,
#             'b': 4,
#             'wettability' : 'Oil-wet'}

# Define relative permeability functions
rel_perm = Relative_perm(kr_infos)

# Define Sw range
Sw = np.linspace(kr_infos['Swc'],1-kr_infos['Sor'],2000)

# Builds and show relative permeability vector values
rel_perm.set_krw(Sw)
rel_perm.set_kro(Sw)
# rel_perm.plot_perm(Sw)

# Definig fractional flow
fw_infos = {'M' : 200.0}  # M = mu_displaced / mu_injected
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

# plt.plot(Sol_vD,Sol_Sw)
# plt.ylim([-0.1,1.1])
# plt.xlim([-0.1,6.0])
# plt.grid(True)
# plt.show()











# Recovery Calculations

Sw1 = np.linspace(Sws, 1-kr_infos['Sor'], 100)
Sw1 = Sw1[1:len(Sw1)-1]

NpD = np.zeros_like(Sw1)
tD1 = np.zeros_like(Sw1)

for i in range(len(Sw1)):
    Sw1_idx = np.argmin(np.abs(Sw - Sw1[i]))
    # print(Sw[Sw1_idx], Sw1[i])
    fw1 = frac_flow.get_fw()[Sw1_idx]
    # print(Sw1[i], fw1)
    tD1[i] = 1 / dfw_dSw[Sw1_idx]
    # print(tD1[i])

    Sw_bar = Sw1[i] + tD1[i]*(1 - fw1)

    NpD[i] = Sw_bar - kr_infos['Swc']

tDmax = 3.0

id_tDmax = np.argmin(np.abs(tDmax - tD1)) + 1

tD = np.concatenate((np.array([0.0, 1/vsD]), tD1[:id_tDmax])) 
NpD = np.concatenate((np.array([0.0, 1/vsD]), NpD[:id_tDmax])) 

plt.plot(tD,NpD,c='b',label='Analítico ' + r'$\mu_{disp.}/\mu_{injec.} = $' + str(fw_infos['M']))
plt.ylim([-0.1, 1.1])
plt.grid(True)
plt.title(kr_infos['wettability'])
plt.ylabel(r"$N_{p_D}$")
plt.xlabel(r"$t_D$")
plt.legend()
plt.show()