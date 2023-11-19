import numpy as np
import matplotlib.pyplot as plt

class Recovery_calc:

    def __init__(self, kr_infos, frac_flow, bl_solution):

        self.kr_infos = kr_infos
        self.bl_sol = bl_solution
        self.fw = frac_flow

    def do_recovery(self):
        Sw1 = np.linspace(self.bl_sol.get_Sws(), 1-self.kr_infos['Sor'], 100)
        Sw1 = Sw1[1:len(Sw1)-1]

        NpD = np.zeros_like(Sw1)
        tD1 = np.zeros_like(Sw1)

        for i in range(len(Sw1)):
            Sw1_idx = np.argmin(np.abs(self.bl_sol.get_Sw() - Sw1[i]))
            # print(Sw[Sw1_idx], Sw1[i])
            fw1 = self.fw.get_fw()[Sw1_idx]
            # print(Sw1[i], fw1)
            tD1[i] = 1 / self.fw.get_dfw_dSw()[Sw1_idx]
            # print(tD1[i])

            Sw_bar = Sw1[i] + tD1[i]*(1 - fw1)

            NpD[i] = Sw_bar - self.kr_infos['Swc']

        self.tD = np.concatenate((np.array([0.0, 1/self.bl_sol.get_vsD()]), tD1)) 
        self.NpD = np.concatenate((np.array([0.0, 1/self.bl_sol.get_vsD()]), NpD)) 

    def show_curve(self, tDmax = 3.0, save_pdf = False):
        # tDmax : max pore volume of injected water to show (nondimensional time)

        id_tDmax = np.argmin(np.abs(tDmax - self.tD))

        plt.plot(self.tD[:id_tDmax], self.NpD[:id_tDmax],c='b',label='Analítico ' + r'$\mu_{disp.}/\mu_{injec.} = $' + str(self.kr_infos['M']))
        plt.ylim([-0.1, 1.1])
        plt.grid(True)
        plt.title('Recovery calculation \n wettability condition: ' + self.kr_infos['wettability'])
        plt.ylabel(r"$N_{p_D}$")
        plt.xlabel(r"$t_D$")
        plt.legend()
        if save_pdf: plt.savefig('recovery_curve_' + self.kr_infos['wettability'] + '.pdf', dpi=300, bbox_inches='tight')
        plt.show()