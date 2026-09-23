# -*- coding: utf-8 -*-
"""
dogma.diagrama
===============
Genera una imagen (PNG) que resume visualmente las tres etapas del
dogma central para una simulación concreta: horquilla de replicación
(hebra líder + fragmentos de Okazaki), transcripción y traducción.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from .simulador import ResultadoReplicacion, ResultadoTranscripcion, ResultadoTraduccion

COLOR_ADN = "#2b6cb0"
COLOR_ADN_NUEVO = "#63b3ed"
COLOR_ARN = "#c05621"
COLOR_PROT = "#2f855a"
COLOR_FONDO_PANEL = "#f7fafc"


def _recortar(seq: str, maximo: int = 34) -> str:
    if len(seq) <= maximo:
        return seq
    mitad = (maximo - 3) // 2
    return seq[:mitad] + "..." + seq[-mitad:]


def _panel(ax, titulo: str):
    ax.set_title(titulo, fontsize=13, fontweight="bold", loc="left", pad=10)
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    ax.add_patch(
        FancyBboxPatch(
            (0.05, 0.05), 9.9, 9.6,
            boxstyle="round,pad=0.02,rounding_size=0.15",
            linewidth=1, edgecolor="#cbd5e0", facecolor=COLOR_FONDO_PANEL,
        )
    )


def generar_diagrama(
    replicacion: ResultadoReplicacion,
    transcripcion: ResultadoTranscripcion,
    traduccion: ResultadoTraduccion,
    ruta_salida: str,
) -> str:
    fig, axes = plt.subplots(3, 1, figsize=(11, 13))
    fig.suptitle(
        "Simulación del Dogma Central de la Biología Molecular",
        fontsize=16, fontweight="bold", y=0.995,
    )

    # ------------------------------------------------------------------
    # Panel 1: Replicación
    # ------------------------------------------------------------------
    ax = axes[0]
    _panel(ax, "1. Replicación del ADN")

    ax.text(0.4, 8.6, "Hebra A (molde hebra líder)", fontsize=9, color="#1a202c")
    ax.plot([0.4, 9.6], [8.1, 8.1], color=COLOR_ADN, lw=3)
    ax.text(0.4, 7.55, "Hebra líder nueva (síntesis continua)", fontsize=9, color="#1a202c")
    ax.plot([0.4, 9.6], [7.1, 7.1], color=COLOR_ADN_NUEVO, lw=3)
    ax.add_patch(FancyArrowPatch((0.4, 7.1), (2.0, 7.1), arrowstyle="-|>",
                                  mutation_scale=14, color=COLOR_ADN_NUEVO))

    ax.text(0.4, 6.3, "Hebra B (molde hebra retrasada)", fontsize=9, color="#1a202c")
    ax.plot([0.4, 9.6], [5.85, 5.85], color=COLOR_ADN, lw=3)

    n_frag = max(1, len(replicacion.fragmentos_okazaki))
    ancho = 9.2 / n_frag
    ax.text(0.4, 5.35, "Fragmentos de Okazaki (síntesis discontinua)", fontsize=9, color="#1a202c")
    for i in range(n_frag):
        x0 = 0.4 + i * ancho
        ax.plot([x0, x0 + ancho * 0.85], [4.9, 4.9], color="#dd6b20", lw=4,
                 solid_capstyle="butt")
        ax.text(x0 + ancho * 0.42, 4.55, f"{i + 1}", fontsize=7, ha="center")
    ax.text(0.4, 4.0, "→ unidos después por la ADN ligasa", fontsize=8,
            color="#4a5568", style="italic")

    ax.text(0.4, 3.1, "Helicasa abre la doble hélice  |  Primasa coloca cebadores"
            "  |  ADN polimerasa III sintetiza  |  Ligasa sella los fragmentos",
            fontsize=8, color="#4a5568", wrap=True)
    ax.text(0.4, 2.3, f"Nº de fragmentos de Okazaki simulados: {len(replicacion.fragmentos_okazaki)}",
            fontsize=8, color="#4a5568")
    ax.text(0.4, 1.5,
            "Resultado: 2 dúplex de ADN hijos, cada uno con una hebra "
            "parental y una hebra nueva (semiconservativa).",
            fontsize=8, color="#4a5568")

    # ------------------------------------------------------------------
    # Panel 2: Transcripción
    # ------------------------------------------------------------------
    ax = axes[1]
    _panel(ax, "2. Transcripción (ADN → ARNm)")
    ax.text(0.4, 8.4, "Hebra molde de ADN (3'→5')", fontsize=9)
    ax.plot([0.4, 9.6], [7.9, 7.9], color=COLOR_ADN, lw=3)
    ax.text(0.4, 7.35, _recortar(transcripcion.hebra_molde[::-1]), fontsize=8, family="monospace")

    ax.add_patch(FancyArrowPatch((5.0, 7.0), (5.0, 5.9), arrowstyle="-|>",
                                  mutation_scale=16, color=COLOR_ARN, lw=2))
    ax.text(5.2, 6.4, "ARN polimerasa", fontsize=8, color=COLOR_ARN, style="italic")

    ax.text(0.4, 5.5, "ARNm sintetizado (5'→3')", fontsize=9)
    ax.plot([0.4, 9.6], [5.0, 5.0], color=COLOR_ARN, lw=3)
    ax.text(0.4, 4.45, _recortar(transcripcion.arnm), fontsize=8, family="monospace")

    ax.text(0.4, 3.4,
            "Complementariedad molde→ARNm:  A→U   T→A   G→C   C→G",
            fontsize=8.5, color="#4a5568")
    ax.text(0.4, 2.6,
            f"Longitud del ARNm: {len(transcripcion.arnm)} nucleótidos",
            fontsize=8, color="#4a5568")

    # ------------------------------------------------------------------
    # Panel 3: Traducción
    # ------------------------------------------------------------------
    ax = axes[2]
    _panel(ax, "3. Traducción (ARNm → proteína)")

    ax.text(0.4, 8.6, "Lectura de codones sobre el ARNm por el ribosoma:", fontsize=9)
    codones = traduccion.codones if traduccion.codones else []
    max_mostrar = 12
    mostrados = codones[:max_mostrar]
    x = 0.4
    for i, cod in enumerate(mostrados):
        es_stop = (cod == traduccion.codon_paro)
        color = "#c53030" if es_stop else "#2c5282"
        ax.add_patch(FancyBboxPatch((x, 7.1), 0.68, 0.7, boxstyle="round,pad=0.02",
                                     facecolor="white", edgecolor=color, linewidth=1.4))
        ax.text(x + 0.34, 7.45, cod, ha="center", va="center", fontsize=7.5,
                family="monospace", color=color)
        x += 0.72
    if len(codones) > max_mostrar:
        ax.text(x + 0.1, 7.45, "...", fontsize=10, va="center")

    ax.text(0.4, 6.3, "Cadena polipeptídica (ribosoma + ARNt + peptidil-transferasa):", fontsize=9)
    aas = traduccion.proteina
    mostrados_aa = aas[:max_mostrar]
    x = 0.4
    for aa in mostrados_aa:
        ax.add_patch(FancyBboxPatch((x, 4.9), 0.68, 0.7, boxstyle="round,pad=0.02",
                                     facecolor=COLOR_PROT, edgecolor=COLOR_PROT, alpha=0.85))
        ax.text(x + 0.34, 5.25, aa.letra1, ha="center", va="center", fontsize=9,
                color="white", fontweight="bold")
        x += 0.72
        if x > 9.0:
            break

    if aas:
        texto_prot = "-".join(a.nombre3 for a in aas)
        if traduccion.codon_paro:
            texto_prot += f"-STOP({traduccion.codon_paro})"
        ax.text(0.4, 3.7, _recortar(texto_prot, 60), fontsize=8, color="#2f855a", family="monospace")
        ax.text(0.4, 3.0, f"Proteína final: {len(aas)} aminoácidos "
                f"({traduccion.secuencia_1letra})", fontsize=8.5, color="#1a202c")
    else:
        ax.text(0.4, 3.7, "No se generó proteína (sin AUG en marco de lectura válido).",
                fontsize=8.5, color="#c53030")

    plt.tight_layout(rect=[0, 0, 1, 0.975])
    fig.savefig(ruta_salida, dpi=170, facecolor="white")
    plt.close(fig)
    return ruta_salida
