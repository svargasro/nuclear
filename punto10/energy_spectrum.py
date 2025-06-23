import matplotlib.pyplot as plt
import numpy as np
import os
import glob
from scipy import stats

def parse_uncertainty(E_str_array, unc_str_array):
    """
    Interpreta la incertidumbre dada como dígitos finales:
      e.g., E_str="79.5104", unc_str="2" -> unc = 2 * 10^(-4) = 0.0002
      Para unc_str="0" -> asigna 1e-8
    """
    returned_sigma = np.zeros(len(E_str_array))
    unc_int_array  = np.zeros(len(E_str_array))

    for index, unc_str in enumerate(unc_str_array):
        unc_int = int(unc_str)
        unc_int_array[index] = unc_int
        E_str = str(E_str_array[index])

        if unc_int == 0:
            unc_int = 1e-7
            unc_int_array[index] = unc_int
            decimals = 0
        else:
            decimals = len(E_str.split('.')[-1])


        returned_sigma[index] = decimals

    return unc_int_array * 10**(-returned_sigma)

def analyze_rotational_spectrum(niveles, element_name, out_dir):
    """
    Grafica E vs J(J+1) con incertidumbres, ajusta linealmente y reporta
    pendiente, intercepto y R^2 con incertidumbre propagada.
    """
    # Preparar datos
    J = np.array([n['J'] for n in niveles])
    x = J * (J + 1)
    y = np.array([n['E'] for n in niveles])

    print(element_name)
    print("J: ",J)
    print("x: ",x)
    print("y: ",y)
    sigma = np.array([n['unc'] for n in niveles], dtype=float)
    sigma = parse_uncertainty(y,sigma)

    # Ajuste ponderado usando polyfit con pesos
    w = 1 / sigma
    fit, cov = np.polyfit(x, y, 1, w=w, cov=True)
    m, b = fit
    m_err, b_err = np.sqrt(np.diag(cov))

    # Cálculo de R^2 y su incertidumbre aproximada
    y_pred = m * x + b
    ss_res = np.sum((y - y_pred)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    R2 = 1 - ss_res / ss_tot

    # Calcular r y error estándar de r
    r, _ = stats.pearsonr(x, y)
    n = len(x)
    se_r = np.sqrt((1 - r**2) / (n - 2))
    R2_err = 2 * abs(r) * se_r

    # Graficar
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.errorbar(x, y, yerr=sigma, fmt='o', ms=2, label='Datos', capsize=1,  ecolor='black')

    x_line = np.linspace(min(x), max(x), 100)
    ax.plot(x_line, m * x_line + b, 'r--', label='Ajuste lineal')
    ax.set_xlabel(r'$J(J+1)$')
    ax.set_ylabel('Energy (keV)')
    ax.set_title(f'Rotational Spectrum: {element_name}')
    ax.legend()

    # Anotar resultados
    textstr = '\n'.join((
        f'm = {m:.3f} ± {m_err:.3f} keV',
        f'b = {b:.1f} ± {b_err:.1f} keV',
        f'$R^2$ = {R2:.3f} ± {R2_err:.3f}'))
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes,
            fontsize=10, va='top', bbox=dict(boxstyle="round,pad=0.3", alpha=0.2))

    # Guardar
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"{element_name}_rot.pdf")
    fig.savefig(out_path)
    plt.close(fig)
    print(f">> Rotational analysis saved: {out_path}")


def plot_nuclear_levels_from_file(txt_path, out_dir=".", threshold=30):

    with open(txt_path, 'r') as f: #Lectura
        lines = [L.strip() for L in f if L.strip()]

    element_name = lines[0]

    niveles = [
        {
            "J":      int(parts[0]),
            "parity": parts[1],
            "E":      float(parts[2]),
            "unc":    int(parts[3]),
        }
        for parts in (line.split() for line in lines[1:])
    ]


    rot_nuclei = {"Gd-158-64", "U-238-92", "Er-166-68"}
    if element_name in rot_nuclei:
        analyze_rotational_spectrum(niveles, element_name, out_dir)



    #Cálculo de offsets para desplazar etiquetas en caso de tener nivles de energía muy cercanos.
    offsets = [0.0] * len(niveles)
    for i, ni in enumerate(niveles):
        for j, nj in enumerate(niveles):
            if i >= j:
                continue
            dE = ni["E"] - nj["E"]
            if abs(dE) < threshold:

                if dE > 0:
                    offsets[i] += threshold * 0.75
                    offsets[j] -= threshold * 0.7
                else:
                    offsets[j] += threshold * 0.75
                    offsets[i] -= threshold * 0.75

    E_max = max(n["E"] for n in niveles) + 100


    fig, ax = plt.subplots(figsize=(4, 6))

    # Grafica todos los niveles con etiquetas ajustadas
    for idx, n in enumerate(niveles):
        E, unc = n["E"], n["unc"]
        label_J = fr"${n['J']}^{n['parity']}$"
        y_label = E + offsets[idx]
        ax.hlines(E, -0.3, 0.3, lw=1)#, color='tab:orange')
        # Etiqueta J^π y energía
        ax.text(-0.5, y_label, label_J, ha='right', va='center',fontsize=8)
        ax.text( 0.5, y_label, f"{E}({unc})", ha='left',  va='center',fontsize=8)

    # Ajustes de estilo
    for spine in ['top', 'right', 'bottom', 'left']:
        ax.spines[spine].set_visible(False)

    # ax.spines['left'].set_visible(True)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(-1, 1)
    ax.set_ylim(-50, E_max)

    # Labels inferiores
    fig.text(0.22, 0.01, r"$\mathbf{J^{\pi}}$", ha='left',  va='center', fontsize=8, fontweight='bold')
    fig.text(0.90, 0.01, "Energy (keV)", ha='right', va='center', fontsize=8, fontweight='bold')

    plt.tight_layout()

    # Guardar
    out_path = os.path.join(out_dir, f"{element_name}.pdf")
    fig.savefig(out_path)
    plt.close(fig)
    print(f">> Guardado: {out_path}")

def plot_all(folder="./info_txt", out_dir="./spectra_pdfs"):
    os.makedirs(out_dir, exist_ok=True)
    for txt_file in glob.glob(os.path.join(folder, "*.txt")): #Recorre todos los archivos txt en info_txt
        plot_nuclear_levels_from_file(txt_file, out_dir=out_dir)



if __name__ == "__main__":
    plot_all()
