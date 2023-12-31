import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

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

    def do_dimensional_NpD_t(self, dimensional_reservoir):
        necessary_keys = ['L', 'phi', 'qt', 'ti', 'tf', 'Nt', 'Nx', 'A']
        for key in necessary_keys:
            assert key in dimensional_reservoir.keys(), f"'{key}' is not defined."

        L = dimensional_reservoir['L']
        A = dimensional_reservoir['A']
        phi = dimensional_reservoir['phi']
        qt = dimensional_reservoir['qt']
        ti = dimensional_reservoir['ti']
        tf = dimensional_reservoir['tf']
        Nt = dimensional_reservoir['Nt']
        Nx = dimensional_reservoir['Nx']

        t = np.divide(self.tD*phi*L,qt)

        
  
        self.Np = self.NpD*phi*L*A # TODO: Veriticar com Rodrigo se multiplica por L e A também para generalizar dimensões reservatório
        self.t = t
    
    def get_t(self):
        return self.t

    def show_dimensional_NpD_t(self,dimensional_reservoir):
        tmax = np.argmin(np.abs(dimensional_reservoir['tf'] - self.t))
        plt.title('Recovery curve - Analytical')
        plt.xlabel(r"$t$")
        plt.ylabel(r"$N_{p}$")
        plt.plot(self.t[:tmax], self.Np[:tmax], c='b')
        plt.grid(True)
        plt.show()

    def define_integrand(self, bl_solution, frac_flow):
        fw_t = np.zeros(bl_solution.get_grid().shape[0])
        for i in range(len(fw_t)):
            Sw = bl_solution.get_grid()[i][-1]
            fw_t[i] = np.interp(Sw, frac_flow.get_Sw(), frac_flow.get_fw())
        self.fw_t = fw_t

    def integrand(self, t, bl_solution):
        return 1 - np.interp(t, bl_solution.get_t(), self.fw_t)
    
    def production_integration(self, bl_solution, dimensional_reservoir):
        prod = []

        for i in range(len(self.t)):
            Integral = quad(self.integrand, 0, self.t[i], args=(bl_solution), limit = 500, epsabs=1.49e-04, epsrel=1.49e-04)[0]
            Integral*=dimensional_reservoir['qt']*dimensional_reservoir['A']
            prod.append(Integral)

        self.prod = np.array(prod)

    def get_prod(self):
        return self.prod
    
    def show_production_prod_t(self, dimensional_reservoir):
        tmax = np.argmin(np.abs(dimensional_reservoir['tf'] - self.t))

        plt.title('Recovery curve - Numerical Integration')
        plt.xlabel(r"$t$")
        plt.ylabel(r"$N_{p}$")
        plt.plot(self.t[:tmax], self.prod[:tmax])
        plt.grid(True)
        plt.show()