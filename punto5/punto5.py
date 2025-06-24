#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

# Parámetros físicos para interacción α-α
V0 = 500.0         # Profundidad en MeV*fm
alpha = 0.7        # fm^-1
m_alpha = 3727.4   # masa de alpha en MeV/c^2
mu = m_alpha / 2   # masa reducida en MeV/c^2
hbarc = 197.3      # MeV·fm

# Definición de potenciales
def V_yukawa(r):
    return -V0 * np.exp(-alpha * r) / r

def V_cent(r, l):
    return (hbarc**2 * l * (l + 1)) / (2 * mu * r**2)

def V_eff(r, l):
    return V_yukawa(r) + V_cent(r, l)

# Rango de distancias (fm)
r = np.linspace(0.1, 3, 500)

# Graficar para ℓ = 0, 2 y 4
for l in [0, 2, 4]:
    plt.figure()
    plt.plot(r, V_yukawa(r), '--', label=r'$V_N$ - Yukawa')
    plt.plot(r, V_cent(r, l), '-', label=fr'$V_{{cent}}$ (ℓ={l})')
    plt.plot(r, V_eff(r, l),  '-.', label=r'$V_{eff}$')
    # plt.axhline(y=0, color='black', alpha=0.8)
    plt.ylim(-800, 800)
    # plt.title(f'Potenciales efectivos para ℓ = {l}') #Quitar para versión final y evitar repetir título
    plt.xlabel('r (fm)')
    plt.ylabel(r'$V_{eff}$ (MeV)')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"p5_potencial_l_{l}.pdf")
