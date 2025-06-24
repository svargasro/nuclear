import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


hbarc = 197.32          # [MeV·fm]
mp = 938.27             # [MeV/c^2] proton's mass
mn = 939.56             # [MeV/c^2] neutron's mass
mu = mp*mn/(mp+mn)      # [MeV/c^2] reduced mass of the n-p system
k = 1.0                 # [fm^-1] transfered wave number


theta = np.linspace(0, np.pi, 1000)     # different scattering angles

# The collision energy comes inside the wavenumber, that is
E_c = 50                           # [MeV] collision energy
k = np.sqrt(2 * mu * E_c) / hbarc  # [fm^-1] wave number

q = 2*k*np.sin(theta/2)                 # [fm^-1] momentum transference


def angulardist(q, V0, alpha):
    return ((2*mu*V0)/(hbarc**2))**2*(1/(alpha**2 + q**2)**2)   # The expression we found for the E.A.D


# Now, I keep  v0 fixed and vary the values of alpha
V0_fixed = 450                  # [MeV·fm]
alpha_vals = [0.3, 0.5, 0.7, 1.2, 1.8]      # [fm^-1] 


sns.set_style("whitegrid")
sns.set_context("paper")

plt.figure(figsize=(10, 6))
for alpha in alpha_vals:
    ds_do = angulardist(q, V0_fixed, alpha)
    plt.plot(np.degrees(theta), ds_do, label=f"$\\alpha$ = {alpha} [fm$^{{-1}}$]")


plt.xlabel("Ángulo θ [grados]")
plt.ylabel(r"$\frac{d\sigma}{d\Omega}$ [fm$^2$/sr]")
plt.yscale("log")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig('dsdov0fix.pdf')
