import numpy as np
import matplotlib.pyplot as plt

# Reference: 
# Martin J. Blunt - Multiphase Flow in Permeable Media_ A Pore-Scale Perspective-Cambridge University Press (2017)
# Section 9.2

class BuckleyLeverettSolution:
    
    def __init__(self, Sw, frac_flow, kr_infos):
        
        self.kr_infos = kr_infos
        self.Sw = Sw
        self.fw = frac_flow
        

    def construct_solution(self):
        # Compute Sws (shock front saturation) and shock solution
        # TODO: Pensar numa abordagem melhor. Está muito gambiarra!
        diff = np.abs(np.divide(self.fw.get_fw(), self.Sw - self.kr_infos['Swc']+1e-6) - self.fw.get_dfw_dSw()) + 1-self.Sw
        self.Sws_idx = np.argmin(diff)
        self.Sws = self.Sw[self.Sws_idx]

        # Compute showck front velocity vsD
        self.vsD = self.fw.get_fw()[self.Sws_idx] / (self.Sws - self.kr_infos['Swc'])         

        # Building solution at the four scenarios 
        # (constant state1, rarefaction, shock front, constant state2)

        # constant state1 (vD < vD_min)
        # Compute solution constant state at vD <= vD_min
        self.vD_min = self.fw.get_dfw_dSw()[-1]
        vD1 = np.linspace(0,self.vD_min,100)
        Sw1 = 1-self.kr_infos['Sor'] * np.ones_like(vD1)

        # rarefaction/showck front
        Sw2 = self.Sw[self.Sws_idx:][::-1]
        Sw2[-1] = self.Sws
        vD2 = self.fw.get_dfw_dSw()[self.Sws_idx:][::-1]
        vD2[-1] = self.vsD

        # constant state2 (vD > vsD)
        vD3 = np.linspace(self.vsD,6.0,500)
        Sw3 = self.kr_infos['Swc'] * np.ones_like(vD3)

        # Whole solution grouping
        self.Sol_vD = np.concatenate((vD1 , vD2, vD3))
        self.Sol_Sw = np.concatenate((Sw1 , Sw2, Sw3))

    def get_solution(self):
        return [self.Sol_vD,self.Sol_Sw]
    
    def get_Sw(self):
        return self.Sw

    def get_Sw(self):
        return self.Sw

    def get_Sws(self):
        return self.Sws

    def get_vsD(self):
        return self.vsD

    def show_solution(self):
        plt.plot(self.Sol_vD, self.Sol_Sw, c='b',label='Analytical solution ' + r'$\mu_{disp.}/\mu_{injec.} = $' + str(self.kr_infos['M']))
        plt.title('Analytical solution \n wettability condition: ' + self.kr_infos['wettability'])
        plt.ylabel(r"$S_w$")
        plt.xlabel(r"$v_D$")
        plt.ylim([-0.1,1.1])
        plt.xlim([-0.1,6.0])
        plt.grid(True)
        plt.legend()
        plt.show()