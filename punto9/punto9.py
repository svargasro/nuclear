# -*- coding: utf-8 -*-
"""Punto9Edna.ipynb

Automatically generated by Colab.

Original file is located at
    https://colab.research.google.com/drive/1Nna-JSSaVK2oW26VrYOUdRFiHjzBB_yq
"""

import numpy as np
import matplotlib.pyplot as plt

# g factors extracted from Wikipedia
g_factors = {'proton':  {'gl': 1, 'gs': 5.585},'neutron': {'gl': 0, 'gs': -3.826}}

# magnetic moments
def mmoment1(gl, gs, j):
    return (j/(j+1))*(gl*(j+3/2)-gs/2)

def mmoment2(gl, gs, j):
    return gl*(j-1/2)+gs/2

# different values of the angular momentum
j_values = np.array([0.5, 1.5, 2.5, 3.5, 4.5])
j_labels = ['1/2', '3/2', '5/2', '7/2', '9/2']

# experimental data extracted from the nuclei table
# (j, experimental magnetic moment)
experimental_data = {
    'proton': [
        (1/2, -0.088402), (5/2, 3.451), (5/2, 2.8215), (1/2, 1.4597),
        (5/2, 3.5432), (7/2, 2.7421), (1/2, 0.002453), (7/2, 4.353),
        (5/2, 4.27545), (5/2, 3.16013), (5/2, 1.382), (1/2, 0.1592),
        (7/2, 2.402), (9/2, 5.5858), (1/2, -0.088402), (1/2, -0.121),
        (9/2, 6.1535), (7/2, 2.54982), (7/2, 2.7014), (9/2, 4.1103),
        (9/2, 5.69445), (1/2, 0.2092), (5/2, 3.357912), (3/2, 2.38633010),
        (3/2, 2.688648910), (3/2, 3.32564268), (3/2, 0.1734),
        (3/2, 1.89109), (3/2, 1.4397565)
    ],
    'neutron': [
        (3/2, -0.3543), (3/2, -0.3748), (1/2, 0.255), (5/2, -0.662),
        (5/2, -0.5833), (3/2, -0.281435), (7/2, 0.0952), (7/2, 0.0952),
        (1/2, -1.00625), (3/2, -0.9641), (1/2, 0.2034), (9/2, -0.39517),
        (9/2, 0.367), (7/2, 0.51), (5/2, -0.68315), (7/2, -0.383),
        (5/2, 0.595), (5/2, 0.5025), (3/2, -0.4038), (1/2, -0.7338),
        (3/2, 0.2705), (9/2, -1.4037), (9/2, -1.4735), (9/2, -1.197),
        (1/2, 0.610), (1/2, 0.601), (1/2, 0.68), (1/2, 0.60089),
        (5/2, 0.7117), (7/2, -1.588), (5/2, 0.6841), (1/2, 0.5058855),
        (1/2, 0.51), (1/2, 0.5414749), (9/2, 0.96), (7/2, -0.786)
    ]
}


def plotline(particle):
    gl = g_factors[particle]['gl']
    gs = g_factors[particle]['gs']
    m1 = mmoment1(gl, gs, j_values)
    m2 = mmoment2(gl, gs, j_values)

    plt.figure()
    plt.plot(j_values, m1, color='black')
    plt.plot(j_values, m2, color='darkorange')
    plt.xticks(j_values, j_labels)
    plt.xlabel('j')
    plt.ylabel('$\\mu$ [J/T]')

    for x, y in experimental_data[particle]:
        plt.scatter(x, y, color='palevioletred')

    if particle == 'proton':
        plt.text(1.2, 3.5, r'$j=l-\frac{1}{2}$', fontsize=12, ha='right', va='bottom')
        plt.text(1.2, 0.5, r'$j=l+\frac{1}{2}$', fontsize=12, ha='right', va='bottom')
    else:
        plt.text(1.2, 1.25, r'$j=l-\frac{1}{2}$', fontsize=12, ha='right', va='bottom')
        plt.text(1.2, -1.8, r'$j=l+\frac{1}{2}$', fontsize=12, ha='right', va='bottom')

    plt.grid(True)
    plt.legend()
    plt.savefig(f'moment_{particle}.pdf')
    plt.show()


plotline('proton')
plotline('neutron')

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Estilo de graficado
sns.set_style("whitegrid")
sns.set_context("paper")

# Factores giromagnéticos (de Wikipedia)
G_FACTORS = {
    'proton':  {'gl': 1.0, 'gs': 5.585},
    'neutron': {'gl': 0.0, 'gs': -3.826}
}

# Momentos magnéticos teóricos (líneas de Schmidt)
def schmidt_lower(gl, gs, j):
    """Caso j = l - 1/2"""
    return (j / (j + 1)) * (gl * (j + 1.5) - 0.5 * gs)

def schmidt_upper(gl, gs, j):
    """Caso j = l + 1/2"""
    return gl * (j - 0.5) + 0.5 * gs

# Valores de j
j_vals = np.array([0.5, 1.5, 2.5, 3.5, 4.5])
j_ticks = ['1/2', '3/2', '5/2', '7/2', '9/2']

# Datos experimentales
EXPERIMENTAL = {
    'proton': [
        (1/2, -0.088402), (5/2, 3.451), (5/2, 2.8215), (1/2, 1.4597),
        (5/2, 3.5432), (7/2, 2.7421), (1/2, 0.002453), (7/2, 4.353),
        (5/2, 4.27545), (5/2, 3.16013), (5/2, 1.382), (1/2, 0.1592),
        (7/2, 2.402), (9/2, 5.5858), (1/2, -0.088402), (1/2, -0.121),
        (9/2, 6.1535), (7/2, 2.54982), (7/2, 2.7014), (9/2, 4.1103),
        (9/2, 5.69445), (1/2, 0.2092), (5/2, 3.357912), (3/2, 2.38633010),
        (3/2, 2.688648910), (3/2, 3.32564268), (3/2, 0.1734),
        (3/2, 1.89109), (3/2, 1.4397565)
    ],
    'neutron': [
        (3/2, -0.3543), (3/2, -0.3748), (1/2, 0.255), (5/2, -0.662),
        (5/2, -0.5833), (3/2, -0.281435), (7/2, 0.0952), (7/2, 0.0952),
        (1/2, -1.00625), (3/2, -0.9641), (1/2, 0.2034), (9/2, -0.39517),
        (9/2, 0.367), (7/2, 0.51), (5/2, -0.68315), (7/2, -0.383),
        (5/2, 0.595), (5/2, 0.5025), (3/2, -0.4038), (1/2, -0.7338),
        (3/2, 0.2705), (9/2, -1.4037), (9/2, -1.4735), (9/2, -1.197),
        (1/2, 0.610), (1/2, 0.601), (1/2, 0.68), (1/2, 0.60089),
        (5/2, 0.7117), (7/2, -1.588), (5/2, 0.6841), (1/2, 0.5058855),
        (1/2, 0.51), (1/2, 0.5414749), (9/2, 0.96), (7/2, -0.786)
    ]
}

def graficar_momentos(particle_type):
    """Grafica líneas de Schmidt y datos experimentales"""
    gl = G_FACTORS[particle_type]['gl']
    gs = G_FACTORS[particle_type]['gs']

    mu_lower = schmidt_lower(gl, gs, j_vals)
    mu_upper = schmidt_upper(gl, gs, j_vals)

    plt.figure(figsize=(8, 5))
    plt.plot(j_vals, mu_lower, label=r'$j = \ell - \frac{1}{2}$', color='steelblue', linestyle='--')
    plt.plot(j_vals, mu_upper, label=r'$j = \ell + \frac{1}{2}$', color='firebrick', linestyle='-.')

    # Datos experimentales
    j_exp, mu_exp = zip(*EXPERIMENTAL[particle_type])
    plt.scatter(j_exp, mu_exp, color='darkgreen', alpha=0.7, edgecolor='black', label='Datos experimentales')

    # Configuración gráfica
    plt.xticks(j_vals, j_ticks)
    plt.xlabel(r'$j$')
    plt.ylabel(r'$\mu$ [$\mu_N$]')
    #plt.title(f'Momentos magnéticos del {particle_type}')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"momento_{particle_type}.pdf", dpi=300)
    plt.show()

# Ejecutar gráficos
graficar_momentos('proton')
graficar_momentos('neutron')