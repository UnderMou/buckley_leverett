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

    def show_solution(self, save_pdf = False):
        plt.plot(self.Sol_vD, self.Sol_Sw, c='b',label='Analytical solution ' + r'$\mu_{disp.}/\mu_{injec.} = $' + str(self.kr_infos['M']))
        plt.title('Analytical solution \n wettability condition: ' + self.kr_infos['wettability'])
        plt.ylabel(r"$S_w$")
        plt.xlabel(r"$v_D$")
        plt.ylim([-0.1,1.1])
        plt.xlim([-0.1,6.0])
        plt.grid(True)
        plt.legend()
        if save_pdf: plt.savefig('buckleyLeverett_solution_' + self.kr_infos['wettability'] + '.pdf', dpi=300, bbox_inches='tight')
        plt.show()

    def do_dimensional_Sw_x(self, dimensional_reservoir):
        necessary_keys = ['L', 'phi', 'qt', 'ti', 'tf', 'Nt', 'Nx']
        for key in necessary_keys:
            assert key in dimensional_reservoir.keys(), f"'{key}' is not defined."

        L = dimensional_reservoir['L']
        phi = dimensional_reservoir['phi']
        qt = dimensional_reservoir['qt']
        ti = dimensional_reservoir['ti']
        tf = dimensional_reservoir['tf']
        Nt = dimensional_reservoir['Nt']
        Nx = dimensional_reservoir['Nx']

        x = np.linspace(0.0, L , Nx)
        t = np.linspace(ti, tf, Nt)

        grid = np.zeros((len(t), len(x)), dtype=float)

        for i in range(grid.shape[0]):
            
            xD = np.divide(x, L)
            tD = np.divide(qt * t[i], phi * L)
            vD = np.divide(xD, tD)

            grid[i][:] = np.interp(vD, self.Sol_vD, self.Sol_Sw)
            
        self.x = x
        self.t = t
        self.grid = grid  

    def get_grid(self):
        return self.grid
    
    def get_t(self):
        return self.t

    def show_dimensional_Sw_x(self):
        # Create the plot
        plt.ion()  # Turn on interactive mode for live updating
        fig, ax = plt.subplots(figsize=(20,5))
        ylimit = [-0.1, 1.1]
        ax.set_ylim(ylimit)
        ax.grid(True)
        ax.set_xlabel('x')
        ax.set_ylabel('Sw')
        line, = ax.plot(self.x, self.grid[0][:])

        text_element = ax.text(0.825, 1.05, '', fontsize=10, ha='left', va='center', color='k')


        for i in range(self.grid.shape[0]):
            new_text = r'$t_D$: ' + f'{self.t[i]:.4f}'
            text_element.set_text(new_text)

            line.set_ydata(self.grid[i][:])   # Update the y-data of the line
            plt.draw()          # Redraw the plot
            plt.pause(0.01)     # Add a small delay to control the update rate
            
        plt.ioff()  # Turn off interactive mode when done
        plt.show()  # Display the final plot

