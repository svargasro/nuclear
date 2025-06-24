import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


A = 43                      # mass number

amu_MeV = 931.5             # 1 amu = 931.5 MeV

# parameters converted from MeV to amu
a_v = 15.753/amu_MeV        # e.g: 15.753 MeV * uma / 931.5 MeV    
a_s = 17.804/amu_MeV
a_c = 0.7103/amu_MeV
a_sym = 23.69/amu_MeV
a_p = 33.6/amu_MeV


mp = 938.27/amu_MeV   # proton's mass converted from MeV to amu
mn = 939.56/amu_MeV   # neutron's mass converted from MeV to amu


def pairing(A, Z):                      # pairing term for the formula
    N = A - Z                           
    if A % 2 == 1:                      # if A is odd, delta=0
        return 0                        
    elif Z % 2 == 0 and N % 2 == 0:     # A even: even-even
        return +a_p / A**(3/4)      
    elif Z % 2 == 1 and N % 2 == 1:     # A even: odd-odd
        return -a_p / A**(3/4)
    

def mass_excess(Z, A):                              
    alpha=mn-a_v+a_s/(A**(1/3))+a_sym               
    beta = (mn-mp)+4*a_sym+a_c/(A**(1/3))
    gamma=4*a_sym/(A)+a_c/(A**(1/3))
    M = alpha*A-beta*Z+gamma*(Z**2)-pairing(A,Z)    # mass formula indicated in the slides
    return M-A                                      # mass excess, the atomic number substracted from the mass formula (in amu)


Z_vals = np.linspace(15, 25, 400)       # some values to put into the formula
mass_vals = mass_excess(Z_vals, A)      

#experimental data collected from https://www-nds.iaea.org/, converted from keV to MeV, and then to amu
Z_exp = np.array([17, 18, 19, 20, 21, 22, 23])
mass_exp_amu = np.array([-24.16, -32.01, -36.57, -38.40, -36.18, -29.31, -17.92]) / amu_MeV


sns.set_style("whitegrid")
sns.set_context("paper")


plt.figure(figsize=(8, 5))
plt.plot(Z_vals, mass_vals)
plt.plot(Z_exp, mass_exp_amu, 'ko', zorder=5)


for i in range(len(Z_exp) - 1):                 # to plot the transitions
    z0, z1 = Z_exp[i], Z_exp[i+1]
    y0, y1 = mass_exp_amu[i], mass_exp_amu[i+1]
    if abs(z1 - z0) == 1:
        if y1 < y0:
            #Bminus
            plt.annotate('', xy=(z1, y1), xytext=(z0, y0),
                         arrowprops=dict(arrowstyle="->", color='black'))
        elif y0 < y1:
            #Bplus
            plt.annotate('', xy=(z0, y0), xytext=(z1, y1),
                         arrowprops=dict(arrowstyle="->", color='black'))

transitions = {                             # these are indicated in the IAEA database
    17: 'β⁻', 18: 'β⁻', 19: 'β⁻',           
    21: 'ec', 22: 'ec', 23: 'ec'
}
for z, label in transitions.items():
    idx = np.where(Z_exp == z)[0]
    if len(idx) > 0:
        y = mass_exp_amu[idx[0]]
        plt.text(z, y + 0.0025, label, fontsize=12, ha='center')

plt.xticks(np.arange(15, 26, 1)) 
plt.xlabel("Z")
plt.ylabel("Exceso de Masa [amu]")
plt.grid(True)
plt.tight_layout()
plt.savefig('A43.pdf')
#plt.show()
