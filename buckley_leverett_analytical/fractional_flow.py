import numpy as np
import matplotlib.pyplot as plt

class Fractional_flow:

    def __init__(self, dict_infos):
        necessary_keys = ['M']
        for key in necessary_keys:
            assert key in dict_infos.keys(), f"'{key}' is not defined."
        
        self.M = dict_infos['M']



    def set_fw(self, Sw, krw, kro):
        np.seterr(all='ignore')
        # Fractional flow ignoring gravity / capillarity effects
        self.fw = np.divide(1, 1 + (1/self.M)*np.divide(kro, krw))
        np.seterr(all='warn')

    def get_fw(self):
        return self.fw

    def plot_fw(self, Sw):
        plt.plot(Sw,self.fw,label='fw')
        plt.grid(True)
        plt.legend()
        plt.ylim([-0.1,1.1])
        plt.xlim([-0.1,1.1])
        plt.show()