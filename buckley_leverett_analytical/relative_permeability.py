import numpy as np
import matplotlib.pyplot as plt

# Reference: 
# Martin J. Blunt - Multiphase Flow in Permeable Media_ A Pore-Scale Perspective-Cambridge University Press (2017)
# Section 7.4.2

class Relative_perm:
    def __init__(self,dict_infos):
        necessary_keys = ['Swc', 'Sor', 'krw_max', 'kro_max', 'a', 'b']
        for key in necessary_keys:
            assert key in dict_infos.keys(), f"'{key}' is not defined."
        
        self.Swc = dict_infos['Swc']
        self.Sor = dict_infos['Sor']
        self.krw_max = dict_infos['krw_max']
        self.kro_max = dict_infos['kro_max']
        self.a = dict_infos['a']
        self.b = dict_infos['b']


    def Se_func(self,Sw):    # Effective water saturation function
        return np.divide(Sw - self.Swc, 1-self.Sor - self.Swc)

    def set_krw(self,Sw):
        self.krw = self.krw_max*np.power(self.Se_func(Sw),self.a)

    def set_kro(self,Sw):
        self.kro = self.kro_max*np.power(1-self.Se_func(Sw),self.b)

    def plot_perm(self, Sw):
        plt.plot(Sw,self.krw,label='krw')
        plt.plot(Sw,self.kro,label='kro')
        plt.grid(True)
        plt.legend()
        plt.ylim([-0.1,1.1])
        plt.xlim([-0.1,1.1])
        plt.show()
