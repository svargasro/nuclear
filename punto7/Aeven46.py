import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


A = 46                      # mass number

amu_MeV = 931.5             # 1 amu = 931.5 MeV

# parameters converted from MeV to amu
a_v = 15.753/amu_MeV        # e.g: 15.753 MeV * uma / 931.5 MeV    
a_s = 17.804/amu_MeV
a_c = 0.7103/amu_MeV
a_sym = 23.69/amu_MeV
a_p = 33.6/amu_MeV


mp = 938.27/amu_MeV   # proton's mass converted from MeV to amu
mn = 939.56/amu_MeV   # neutron's mass converted from MeV to amu
    

def mass_excess(Z, A, parity):              # the mass formula indicated in the slides                
    alpha=mn-a_v+a_s/(A**(1/3))+a_sym               
    beta = (mn-mp)+4*a_sym+a_c/(A**(1/3))
    gamma=4*a_sym/(A)+a_c/(A**(1/3))

    if parity == 'even-even':               # for A even: even-even pairing
        delta = +a_p / A**(3/4)             
    elif parity == 'odd-odd':               # for A even: odd-odd pairing
        delta = -a_p / A**(3/4)
    else:
        delta = 0                           # for A odd

    M = alpha*A-beta*Z+gamma*(Z**2)-delta           
    return M-A                                      # mass excess (in amu)


Z_vals = np.linspace(16, 26, 400)       # some values to put into the formula
mass_vals_even = [mass_excess(Z, A, 'even-even') for Z in Z_vals] # even-even
mass_vals_odd = [mass_excess(Z, A, 'odd-odd') for Z in Z_vals]    # odd-odd

#experimental data collected from https://www-nds.iaea.org/, converted from keV to MeV and then to amu
Z_exp = np.array([18, 19, 20, 21, 22, 23, 24])
mass_exp_amu = np.array([-29.77, -35.41, -43.13, -41.76, -44.12, -37.07, -29.47]) / amu_MeV


sns.set_style("whitegrid")
sns.set_context("paper")


plt.figure(figsize=(8, 5))
plt.plot(Z_vals, mass_vals_even, label='par-par', linestyle='-')
plt.plot(Z_vals, mass_vals_odd, label='impar-impar', linestyle='-')
plt.plot(Z_exp, mass_exp_amu,'ko', zorder=5, label='exp data')


for i in range(len(Z_exp) - 1):                 # to plot the transitions
    z0, z1 = Z_exp[i], Z_exp[i+1]
    y0, y1 = mass_exp_amu[i], mass_exp_amu[i+1]
    if abs(z1 - z0) == 1:
        if (z0 == 20 and z1 == 21):         
            continue                            # I skip the arrow because Sc doesnt decay to Ca
        if y1 < y0:
            #Bminus
            plt.annotate('', xy=(z1, y1), xytext=(z0, y0),
                         arrowprops=dict(arrowstyle="->", color='black'))
        elif y0 < y1:
            #Bplus
            plt.annotate('', xy=(z0, y0), xytext=(z1, y1),
                         arrowprops=dict(arrowstyle="->", color='black'))

beta_transitions = {
    18: 'β⁻', 19: 'β⁻',
    21: 'β⁻', 23: 'ec', 24: 'ec'
}
for z, label in beta_transitions.items():
    idx = np.where(Z_exp == z)[0]
    if len(idx) > 0:
        y = mass_exp_amu[idx[0]]
        plt.text(z, y + 0.0025, label, fontsize=12, color='black', ha='center')

plt.xticks(np.arange(16, 26, 1)) 
plt.xlabel("Z")
plt.ylabel("Exceso de Masa [amu]")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('A46.pdf')
#plt.show()
