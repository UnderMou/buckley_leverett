import numpy as np
import matplotlib.pyplot as plt

class Fractional_flow:

    def __init__(self, Sw, rel_perm, dict_infos):
        necessary_keys = ['M']
        for key in necessary_keys:
            assert key in dict_infos.keys(), f"'{key}' is not defined."
        
        self.kr_infos = dict_infos
        self.Sw = Sw
        self.rel_perm = rel_perm
        self.M = dict_infos['M']

    def set_fw(self):
        krw, kro = self.rel_perm.get_krw(), self.rel_perm.get_kro()
        np.seterr(all='ignore')

        # Fractional flow ignoring gravity / capillarity effects
        self.fw = np.divide(1, 1 + (1/self.M)*np.divide(kro, krw))

        np.seterr(all='warn')
        self.dfw_dSw = np.gradient(self.fw, self.Sw)

    def get_fw(self):
        return self.fw
    
    def get_dfw_dSw(self):
        return self.dfw_dSw

    def get_Sw(self):
        return self.Sw

    def plot_fw(self, save_pdf = False):
        plt.title('Fractional flow')
        plt.xlabel(r"$S_w$")
        plt.ylabel(r"$f_w$")
        plt.plot(self.Sw,self.fw,c='k',label='fw')
        plt.grid(True)
        plt.legend()
        plt.ylim([-0.1,1.1])
        plt.xlim([-0.1,1.1])
        if save_pdf: plt.savefig('fractional_flow_' + self.kr_infos['wettability'] + '.pdf', dpi=300, bbox_inches='tight')
        plt.show()

    def eval_fw(self, Sw_eval):
        return np.interp(Sw_eval, self.Sw, self.fw)
    
class Fractional_flow_grav(Fractional_flow):

    def __init__(self, Sw, rel_perm, dict_infos):
        super().__init__(Sw, rel_perm, dict_infos) 

        self.muw = dict_infos['muw']
        self.rho_w = dict_infos['rho_w']
        self.rho_o = dict_infos['rho_o']
        self.qt = dict_infos['qt']
        self.K = dict_infos['K']
        self.g = dict_infos['g']

    def set_fw(self):
        krw, kro = self.rel_perm.get_krw(), self.rel_perm.get_kro()
        np.seterr(all='ignore')

        # Fractional flow ignoring just capillarity effects
        muw = self.muw
        muo = self.M*muw
        rhow = self.rho_w
        rhoo = self.rho_o
        qt = self.qt
        K = self.K
        g = self.g

        mobw = krw/muw
        mobo = kro/muo
        mobt = mobo + mobw
        self.fw = np.divide(mobw,mobt) + (K/qt)*np.divide(np.multiply(mobw,mobo),mobt)*(rhow-rhoo)*g

        np.seterr(all='warn')
        self.dfw_dSw = np.gradient(self.fw, self.Sw)
