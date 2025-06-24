import matplotlib.pyplot as plt
from matplotlib.patches import Circle

# Desactivar uso de TeX externo para que funcione con mathtext
plt.rcParams['text.usetex'] = False

# Datos de niveles y posiciones verticales con etiquetas en mathtext
levels = [
    (r"$1s_{1/2}$", 0),
    (r"$1p_{3/2}$", 1),
    (r"$1p_{1/2}$", 2),
    (r"$1d_{5/2}$", 3),
    (r"$2s_{1/2}$", 4)
]

data = {
    "14N": {
        "protons": {
            r"$1s_{1/2}$": [+1, -1],
            r"$1p_{3/2}$": [+1, -1, +1, -1],
            r"$1p_{1/2}$": [+1]
        },
        "neutrons": {
            r"$1s_{1/2}$": [+1, -1],
            r"$1p_{3/2}$": [+1, -1, +1, -1],
            r"$1p_{1/2}$": [+1]
        }
    },
    "17O": {
        "protons": {
            r"$1s_{1/2}$": [+1, -1],
            r"$1p_{3/2}$": [+1, -1, +1, -1],
            r"$1p_{1/2}$": [+1, -1]
        },
        "neutrons": {
            r"$1s_{1/2}$": [+1, -1],
            r"$1p_{3/2}$": [+1, -1, +1, -1],
            r"$1p_{1/2}$": [+1, -1],
            r"$1d_{5/2}$": [+1]
        }
    },
    "32S": {
        "protons": {
            r"$1s_{1/2}$": [+1, -1],
            r"$1p_{3/2}$": [+1, -1, +1, -1],
            r"$1p_{1/2}$": [+1, -1],
            r"$1d_{5/2}$": [+1, -1, +1, -1, +1, -1],
            r"$2s_{1/2}$": [+1, -1]
        },
        "neutrons": {
            r"$1s_{1/2}$": [+1, -1],
            r"$1p_{3/2}$": [+1, -1, +1, -1],
            r"$1p_{1/2}$": [+1, -1],
            r"$1d_{5/2}$": [+1, -1, +1, -1, +1, -1],
            r"$2s_{1/2}$": [+1, -1]
        }
    },
    "31P": {
        "protons": {
            r"$1d_{5/2}$": [+1, -1, +1, -1, +1, -1],
            r"$2s_{1/2}$": [+1]
        },
        "neutrons": {
            r"$1s_{1/2}$": [+1, -1],
            r"$1p_{3/2}$": [+1, -1, +1, -1],
            r"$1p_{1/2}$": [+1, -1],
            r"$1d_{5/2}$": [+1, -1, +1, -1, +1, -1],
            r"$2s_{1/2}$": [+1, -1]
        }
    }
}

def draw_and_save_centered(nucleus, side):
    fig, ax = plt.subplots(figsize=(6, 4))
    sym = "π" if side == "protons" else "ν"
    ax.set_title(f"{nucleus} ({sym})", fontsize=14)
    ax.set_xlim(0, 6)
    ax.set_ylim(-0.5, 4.5)
    ax.axis('off')

    # Dibujar niveles con mathtext
    for label, y in levels:
        ax.hlines(y, 0.5, 5.5, color='orange')
        ax.text(5.6, y, label, va='center', fontsize=12)
    ax.text(0.3, 4.8, sym, fontsize=16)

    # Dibujar ocupaciones centradas
    occ = data[nucleus][side]
    for label, y in levels:
        if label in occ:
            spins = occ[label]
            k = len(spins)
            spacing = 0.6
            xs = [3 + (i - (k-1)/2)*spacing for i in range(k)]
            for x, spin in zip(xs, spins):
                circ = Circle((x, y), 0.2, fill=True, facecolor='gray', edgecolor='black')
                ax.add_patch(circ)
                dy = 0.3 if spin > 0 else -0.3
                ax.arrow(x, y, 0, dy, head_width=0.07, length_includes_head=True, color='black')

    filename = f"{nucleus}_{side}.pdf"
    fig.savefig(f"./{filename}", bbox_inches='tight')
    plt.close(fig)
    return filename

saved = []
for nucleus in data:
    saved.append(draw_and_save_centered(nucleus, "protons"))
    saved.append(draw_and_save_centered(nucleus, "neutrons"))
