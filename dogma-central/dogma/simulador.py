"""Integra las tres etapas del dogma central: ADN → ADN → ARN → proteína."""

from __future__ import annotations

from .replicacion import replicar
from .secuencias import contenido_gc, limpiar_secuencia
from .traduccion import traducir
from .transcripcion import transcribir

LONGITUD_MINIMA = 6
LONGITUD_MAXIMA = 300


def simular(
    secuencia: str,
    tam_fragmento: int = 8,
    tam_cebador: int = 3,
    hebra_molde: str = "inferior",
) -> dict:
    """Ejecuta el flujo completo a partir de la hebra superior (5'→3') de un ADN.

    La transcripción se hace sobre la molécula hija 1 obtenida en la
    replicación, de modo que la información pasa realmente de una etapa a la
    siguiente.
    """
    adn = limpiar_secuencia(secuencia)
    if len(adn) < LONGITUD_MINIMA:
        raise ValueError(f"La secuencia debe tener al menos {LONGITUD_MINIMA} nucleótidos.")
    if len(adn) > LONGITUD_MAXIMA:
        raise ValueError(
            f"La secuencia tiene {len(adn)} nt; el máximo es {LONGITUD_MAXIMA} "
            "para que la visualización sea legible."
        )

    replicacion = replicar(adn, tam_fragmento=tam_fragmento, tam_cebador=tam_cebador)
    hija = replicacion["hijas"][0]["superior"]
    transcripcion = transcribir(hija, hebra_molde=hebra_molde)
    traduccion = traducir(transcripcion["arnm"])

    return {
        "entrada": {"secuencia": adn, "longitud": len(adn), "gc": contenido_gc(adn)},
        "replicacion": replicacion,
        "transcripcion": transcripcion,
        "traduccion": traduccion,
    }
