import matplotlib.pyplot as plt
import numpy as np
import os
import glob
from scipy import stats
import seaborn as sns

def fmt_with_uncertainty(val, err, unit=""):
    """
    Devuelve un string 'v±u[ unit]' donde:
     - u tiene 1 cifra significativa.
     - v está redondeado al mismo decimal que u.
     - Si no hay decimales, no muestra '.0'.
    """
    # Si err es cero muestre sin decimales
    if err == 0:
        return f"{int(val)}±0{unit}"
    # Orden de magnitud de la incertidumbre
    exp = int(np.floor(np.log10(abs(err))))
    # Incertidumbre con 1 cifra significativa
    err_1sig = round(err, -exp)
    # Número de decimales
    dec = max(-exp, 0)
    # Elegir formato según dec
    if dec > 0:
        fmt = f"{{val:.{dec}f}}±{{err:.{dec}f}}{{unit}}"
    else:
        fmt = f"{{val:.0f}}±{{err:.0f}}{{unit}}"
    return fmt.format(val=val, err=err_1sig, unit=unit)

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

    J = np.array([n['J'] for n in niveles])
    x = J * (J + 1)
    y = np.array([n['E'] for n in niveles])

    # print(element_name)
    # print("J: ",J)
    # print("x: ",x)
    # print("y: ",y)
    sigma = np.array([n['unc'] for n in niveles], dtype=float)
    sigma = parse_uncertainty(y,sigma)

    # Ajuste ponderado usando polyfit con pesos
    w = 1 / sigma
    fit, cov = np.polyfit(x, y, 1, w=w, cov=True)
    m, b = fit
    m_err, b_err = np.sqrt(np.diag(cov))

    # Cálculo de R^2
    y_pred = m * x + b
    ss_res = np.sum((y - y_pred)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    R2 = 1 - ss_res / ss_tot


    r, _ = stats.pearsonr(x, y)
    n = len(x)
    se_r = np.sqrt((1 - r**2) / (n - 2))
    R2_err = 2 * abs(r) * se_r


    sns.set_style("whitegrid")
    sns.set_context("paper")

    # Graficar
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.errorbar(x, y, yerr=sigma, fmt='o', ms=3, capsize=1,  ecolor='black')
    print(element_name)
    print("Máx error: ", np.max(sigma))

    x_line = np.linspace(min(x), max(x), 100)
    ax.plot(x_line, m * x_line + b, '-', label=r'Ajuste lineal: $E=m \cdot J(J+1)+b$ ')
    ax.set_xlabel(r'$J(J+1)$')
    ax.set_ylabel('Energía (keV)')
    # ax.set_title(f'Rotational Spectrum: {element_name}')
    ax.legend()




    # textstr = '\n'.join((
    #     f'm = {m:.3f} ± {m_err:.3f} keV',
    #     f'b = {b:.1f} ± {b_err:.1f} keV',
    #     f'$R^2$ = {R2:.3f} ± {R2_err:.3f}'))
    # ax.text(0.6, 0.25, textstr, transform=ax.transAxes, fontsize=8, va='top')

    # Momento de inercia I = ħ^2/(2m)
    # ħc = 197300 keV·fm
    hbarc = 197300  # keV·fm
    I = (hbarc**2) / (2 * m)           # keV·fm^2/c^2
    I_err = I * (m_err / m)

    m_str  = fmt_with_uncertainty(m,   m_err,   " keV")
    b_str  = fmt_with_uncertainty(b,   b_err,   " keV")
    I_str = fmt_with_uncertainty(I, I_err,"")

    textstr = "\n".join([
        rf"$m = {m_str}$",
        rf"$b = {b_str}$",
        rf"$R^2 = {R2:.3f}$",
        rf"$I = {I_str} \frac{{keV}}{{c^2}}·fm^2$"
    ])
    ax.text(0.55, 0.25, textstr, transform=ax.transAxes,
            fontsize=8, va='top')

    # Guardar
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"rot_{element_name}.pdf")
    fig.savefig(out_path)
    plt.close(fig)
    print(f">> Análisis rotational guardado: {out_path}")


def analyze_vibrational_spectrum(niveles, element_name, out_dir):


    datos = sorted(niveles, key=lambda n: n['E'])
    E_vals = np.array([d['E'] for d in datos])
    unc_ints = np.array([d['unc'] for d in datos])
    sigmas = parse_uncertainty(E_vals, unc_ints)


    # E0, E1 y delta
    E0, E1 = E_vals[:2]
    sigma0, sigma1 = sigmas[:2]
    delta = (E1 - E0)*0.2

    # Organización para promedios
    clusters = []
    errors = []
    if len(E_vals) > 2:
        current = [E_vals[2]]
        current_sig = [sigmas[2]]
        for prev, curr, sig in zip(E_vals[2:], E_vals[3:], sigmas[3:]):
            if abs(curr - prev) < delta: #Si están muy cerca, se añaden al mismo cluster.
                current.append(curr)
                current_sig.append(sig)
            else:
                clusters.append(current)
                errors.append(np.sqrt(np.sum(np.array(current_sig)**2)) / len(current_sig)) #Error con desviación estándar sobre N
                current = [curr]
                current_sig = [sig]
        clusters.append(current)
        errors.append(np.sqrt(np.sum(np.array(current_sig)**2)) / len(current_sig))


    n_vals = [0, 1] + list(range(2, 2 + len(clusters)))
    En = [E0, E1] + [np.mean(c) for c in clusters] #Promedios de cada cluster
    sigma_n = [sigma0, sigma1] + errors

    #Se ignora el valor de 0,0.
    En = En[1:]
    n_vals = n_vals[1:]
    sigma_n = sigma_n[1:]

    # Ajuste ponderado En vs n
    x = np.array(n_vals)
    y = np.array(En)
    sigma_n = np.array(sigma_n)
    w = 1 / sigma_n
    fit, cov = np.polyfit(x, y, 1, w=w, cov=True)
    m, b = fit
    m_err, b_err = np.sqrt(np.diag(cov))

    # Cálculo R2 y su incertidumbre
    y_pred = m * x + b
    ss_res = np.sum((y - y_pred)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    R2 = 1 - ss_res / ss_tot



    # print(element_name)
    # print("Promedidos: ", y)
    # print("Total: ", E_vals)


    # Graficar
    fig, ax = plt.subplots(figsize=(5,4))
    ax.errorbar(n_vals, En, yerr=sigma_n, fmt='o', ms=2, capsize=1, ecolor="black")
    print(element_name)
    print("Máx error: ", np.max(sigma_n))
    x_line = np.linspace(0, n_vals[-1], 100)
    ax.plot(x_line, m*x_line + b, '-', label=r'Ajuste lineal: $E_n = m \cdot n+b$')
    ax.set_xticks(n_vals)
    ax.set_xlabel('n')
    ax.set_ylabel(r'$E_n$ (KeV)')
    # ax.set_title(f'Modelo vibracional: {element_name}')
    ax.legend()

    # Preparar texto con resultados
    # textstr = '\n'.join([
    #     rf'$m = {m:.1f}\pm{m_err:.1f}\ \mathrm{{keV}}$,
    #     rf'$b = {b:.4f}\pm{b_err:.4f}\ \mathrm{{keV}}$',
    #     rf'$R^2 = {R2:.3f}$',
    #     # rf'$\omega = \dfrac{{m}}{{\hbar}} \approx {omega:.3e}\ \mathrm{{s}}^{{-1}}$'
    # ])


    m_str  = fmt_with_uncertainty(m,   m_err,   " keV")
    b_str  = fmt_with_uncertainty(b,   b_err,   " keV")
    omega_str  = fmt_with_uncertainty(m,   m_err,   "")

    textstr = "\n".join([
        rf"$m = {m_str}$",
        rf"$R^2 = {R2:.3f}$",
        rf"$\omega = {omega_str}$  $\frac{{keV}}{{\hbar}}$"
    ])
    ax.text(0.6, 0.25, textstr, transform=ax.transAxes,
            fontsize=8, va='top')

    # Guardar
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f"vib_{element_name}.pdf")
    plt.savefig(path)
    plt.close()
    print(f">> Análisis vibracional guardado: {path}")


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

    vib_nuclei = {"Cd-110-48", "Pd-106-46", "Zr-90-40"}
    if element_name in vib_nuclei:
        analyze_vibrational_spectrum(niveles, element_name, out_dir)

    #Cálculo de offsets para desplazar etiquetas en caso de tener nivles de energía muy cercanos.
    offsets = [0.0] * len(niveles)
    for i, ni in enumerate(niveles):
        for j, nj in enumerate(niveles):
            if i >= j:
                continue
            dE = ni["E"] - nj["E"]
            if abs(dE) < threshold:

                if dE > 0:
                    offsets[i] += threshold * 0.65
                    offsets[j] -= threshold * 0.65
                else:
                    offsets[j] += threshold * 0.65
                    offsets[i] -= threshold * 0.65

    E_max = max(n["E"] for n in niveles) + 100



    sns.set_style("whitegrid")
    sns.set_context("paper")

    fig, ax = plt.subplots(figsize=(4, 6))

    # Grafica de todos los niveles con etiquetas ajustadas
    for idx, n in enumerate(niveles):
        E, unc = n["E"], n["unc"]
        label_J = fr"${n['J']}^{n['parity']}$"
        y_label = E + offsets[idx]
        ax.hlines(E, -0.3, 0.3, lw=0.5)#, color='tab:orange')
        # Etiqueta J^π y energía
        ax.text(-0.5, y_label, label_J, ha='right', va='center',fontsize=8)
        ax.text( 0.5, y_label, f"{E}({unc})", ha='left',  va='center',fontsize=8)

    # Ajustes de estilo
    for spine in ['top', 'right', 'bottom', 'left']:
        ax.spines[spine].set_visible(True)

    # ax.spines['left'].set_visible(True)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(-1, 1)
    ax.set_ylim(-50, E_max)

    # Labels inferiores
    fig.text(0.22, 0.01, r"$\mathbf{J^{\pi}}$", ha='left',  va='center', fontsize=8, fontweight='bold')
    fig.text(0.90, 0.01, "Energía (keV)", ha='right', va='center', fontsize=8, fontweight='bold')

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
