import matplotlib.pyplot as plt
import os
import glob

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
