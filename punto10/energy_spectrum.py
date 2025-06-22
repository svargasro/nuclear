import matplotlib.pyplot as plt
import os
import glob

def plot_nuclear_levels_from_file(txt_path, out_dir="."):
    # Leer líneas no vacías
    with open(txt_path, 'r') as f:
        lines = [L.strip() for L in f if L.strip()]


    element_name = lines[0]     # Primera línea: nombre del núcleo
    # Resto: parseo con comprensión de listas
    niveles = [
        {
            "J":       int(parts[0]),
            "parity":  parts[1],
            "E":       float(parts[2]),
            "unc":     int(parts[3]),
        }
        for parts in (line.split() for line in lines[1:])
    ]

    # Precalcular máximos
    E_max = max(n["E"] for n in niveles) + 100

    # Crear figura
    fig, ax = plt.subplots(figsize=(4, 8))
    # ax.axvline(0, color='black', lw=1)

    # Graficar todos los niveles
    for n in niveles:
        E, unc = n["E"], n["unc"]
        label_J = fr"${n['J']}^{n['parity']}$"
        ax.hlines(E, -0.3, 0.3, lw=1, color='tab:orange')
        ax.text(-0.5, E, label_J, ha='right', va='center')
        ax.text( 0.5, E, f"{E:.1f}({unc})", ha='left',  va='center')

    # Ajustes de estilo: solo la línea de referencia y nada más
    for spine in ['top', 'right', 'bottom', 'left']:
        ax.spines[spine].set_visible(False)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(-1, 1)
    ax.set_ylim(-50, E_max)

    # Labels inferiores
    fig.text(0.22, 0.02, r"$J^{\pi}$", ha='left',  va='center', fontsize=10)
    fig.text(0.90, 0.02, "Energy (keV)", ha='right', va='center', fontsize=10)

    plt.tight_layout()

    # Guardar
    out_path = os.path.join(out_dir, f"{element_name}.pdf")
    fig.savefig(out_path)
    plt.close(fig)
    print(f">> Guardado: {out_path}")

def plot_all(folder="./info_txt", out_dir="."):
    os.makedirs(out_dir, exist_ok=True)
    for txt in glob.glob(os.path.join(folder, "*.txt")):
        plot_nuclear_levels_from_file(txt, out_dir=out_dir)

if __name__ == "__main__":
    plot_all(folder="./info_txt", out_dir="./spectra_pdfs")
