import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

hbarc = 197.32          # [MeV·fm]
mp = 938.27             # [MeV/c^2] proton's mass
mn = 939.56             # [MeV/c^2] neutron's mass
mu = mp*mn/(mp+mn)      # [MeV/c^2] reduced mass of the n-p system

# The collision energy comes inside the wavenumber, that is
E_c = 50                           # [MeV] collision energy
k = np.sqrt(2 * mu * E_c) / hbarc  # [fm^-1] wave number


theta = np.linspace(0, np.pi, 1000)     # different scattering angles
q = 2*k*np.sin(theta/2)                 # [fm^-1] momentum transference


def angulardist(q, V0, alpha):
    return ((2*mu*V0)/(hbarc**2))**2*(1/(alpha**2 + q**2)**2)   # The expression we found for the E.A.D


alpha_fixed = 0.7                 # [fm^-1]
V0_vals = [20, 50, 80, 100, 200, 300, 600]  # [MeV·fm]


sns.set_style("whitegrid")
sns.set_context("paper")


plt.figure(figsize=(10, 6))
for V0 in V0_vals:
    ds_do = angulardist(q, V0, alpha_fixed)
    plt.plot(np.degrees(theta), ds_do, label=f"$V_0$ = {V0} [MeV·fm]")


plt.xlabel("Ángulo θ [grados]")
plt.ylabel(r"$\frac{d\sigma}{d\Omega}$ [fm$^2$/sr]")
plt.yscale('log')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('dsdoafixed.pdf')
