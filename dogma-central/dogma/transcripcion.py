"""Simulación de la transcripción (ADN → ARN).

La ARN polimerasa lee la hebra molde en sentido 3'→5' y sintetiza el ARN
mensajero en sentido 5'→3'. El ARNm resultante es complementario a la hebra
molde e igual a la hebra codificante, con uracilo (U) en lugar de timina (T).

* Si la hebra molde es la **inferior**, la polimerasa avanza de izquierda a
  derecha y el ARNm coincide con la hebra superior.
* Si la hebra molde es la **superior**, avanza de derecha a izquierda y el
  ARNm es la complementaria inversa de la hebra superior.
"""

from __future__ import annotations

import math

from .secuencias import arn_desde_molde, complementar

PASOS_ELONGACION = 18


def transcribir(superior: str, hebra_molde: str = "inferior") -> dict:
    if hebra_molde not in ("inferior", "superior"):
        raise ValueError("La hebra molde debe ser 'inferior' o 'superior'.")

    n = len(superior)
    inferior = complementar(superior)

    if hebra_molde == "inferior":
        molde, codificante = inferior, superior
        arn_alineado = arn_desde_molde(inferior)
        arnm = arn_alineado
        sentido = "derecha"
    else:
        molde, codificante = superior, inferior
        arn_alineado = arn_desde_molde(superior)
        arnm = arn_alineado[::-1]
        sentido = "izquierda"

    fotos = [
        {
            "etapa": "iniciacion",
            "titulo": "Iniciación",
            "texto": (
                "La ARN polimerasa (con los factores de transcripción en eucariotas, o el "
                "factor σ en bacterias) se une al promotor y separa las dos hebras. Se "
                "forma la burbuja de transcripción. Solo la hebra "
                f"{hebra_molde} sirve de molde; la otra, la hebra codificante, tiene la "
                "misma secuencia que tendrá el ARNm (con T en lugar de U). A diferencia "
                "de la ADN polimerasa, la ARN polimerasa no necesita cebador."
            ),
            "transcritos": 0,
        }
    ]

    paso = max(1, math.ceil(n / PASOS_ELONGACION))
    k = 0
    while k < n:
        nuevo = min(n, k + paso)
        tramo = arnm[k:nuevo]
        fotos.append(
            {
                "etapa": "elongacion",
                "titulo": f"Elongación ({nuevo}/{n} nt)",
                "texto": (
                    "La ARN polimerasa lee la hebra molde 3'→5' y añade ribonucleótidos "
                    f"trifosfato al extremo 3' del ARN: 5'…{tramo}-3'. Empareja A con U, T con A, "
                    "C con G y G con C. Detrás de la polimerasa, el ARN se separa del molde "
                    "y la doble hélice se vuelve a cerrar."
                ),
                "transcritos": nuevo,
            }
        )
        k = nuevo

    fotos.append(
        {
            "etapa": "terminacion",
            "titulo": "Terminación",
            "texto": (
                "La ARN polimerasa llega a la señal de terminación, se separa del ADN y "
                f"libera el ARNm de {n} nt. El ADN recupera su doble hélice. En eucariotas, "
                "el pre-ARNm se procesa antes de salir del núcleo: se le añade la caperuza "
                "5' y la cola poli-A, y el corte y empalme elimina los intrones. Este paso "
                "no se simula."
            ),
            "transcritos": n,
        }
    )

    return {
        "hebra_molde": hebra_molde,
        "molde": molde,
        "codificante": codificante,
        "superior": superior,
        "inferior": inferior,
        "arn_alineado": arn_alineado,
        "arnm": arnm,
        "sentido": sentido,
        "fotos": fotos,
    }
