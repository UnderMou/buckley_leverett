import numpy as np
import matplotlib.pyplot as plt

# Reference: 
# Martin J. Blunt - Multiphase Flow in Permeable Media_ A Pore-Scale Perspective-Cambridge University Press (2017)
# Section 7.4.2

class Relative_perm:

    def __init__(self, Sw, dict_infos):
        necessary_keys = ['Swc', 'Sor', 'krw_max', 'kro_max', 'a', 'b']
        for key in necessary_keys:
            assert key in dict_infos.keys(), f"'{key}' is not defined."
        
        self.Swc = dict_infos['Swc']
        self.Sor = dict_infos['Sor']
        self.krw_max = dict_infos['krw_max']
        self.kro_max = dict_infos['kro_max']
        self.a = dict_infos['a']
        self.b = dict_infos['b']
        self.Sw = Sw
        self.case = dict_infos['wettability']



    def Se_func(self):    # Effective water saturation function
        return np.divide(self.Sw - self.Swc, 1-self.Sor - self.Swc)

    def set_krw(self):
        self.krw = self.krw_max*np.power(self.Se_func(),self.a)

    def set_kro(self):
        self.kro = self.kro_max*np.power(1-self.Se_func(),self.b)

    def get_krw(self):
        return self.krw
    
    def get_kro(self):
        return self.kro

    def plot_perm(self, save_pdf = False):
        plt.title('Relative permeability')
        plt.xlabel(r"$S_w$")
        plt.ylabel(r"$k_r$")
        plt.plot(self.Sw,self.krw,c='b',label='krw')
        plt.plot(self.Sw,self.kro,c='k',label='kro')
        plt.grid(True)
        plt.legend()
        plt.ylim([-0.1,1.1])
        plt.xlim([-0.1,1.1])
        if save_pdf: plt.savefig('rel_perm_' + self.case + '.pdf', dpi=300, bbox_inches='tight')
        plt.show()

