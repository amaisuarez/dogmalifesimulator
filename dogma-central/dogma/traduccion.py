"""Simulación de la traducción (ARN → proteína).

El ribosoma recorre el ARNm en sentido 5'→3' hasta el primer codón AUG. Desde
ahí lee tripletes (codones) sin solapamiento. Cada codón se empareja con el
anticodón de un ARNt cargado con su aminoácido, y la cadena crece hasta que
un codón de parada (UAA, UAG, UGA) entra en el sitio A.
"""

from __future__ import annotations

from itertools import product

from .secuencias import complementar_arn

# Código genético estándar, en el orden clásico de la tabla (U, C, A, G).
_BASES = "UCAG"
_AMINOACIDOS_1 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"

AMINOACIDOS = {
    "A": ("Ala", "Alanina"),
    "R": ("Arg", "Arginina"),
    "N": ("Asn", "Asparagina"),
    "D": ("Asp", "Ácido aspártico"),
    "C": ("Cys", "Cisteína"),
    "Q": ("Gln", "Glutamina"),
    "E": ("Glu", "Ácido glutámico"),
    "G": ("Gly", "Glicina"),
    "H": ("His", "Histidina"),
    "I": ("Ile", "Isoleucina"),
    "L": ("Leu", "Leucina"),
    "K": ("Lys", "Lisina"),
    "M": ("Met", "Metionina"),
    "F": ("Phe", "Fenilalanina"),
    "P": ("Pro", "Prolina"),
    "S": ("Ser", "Serina"),
    "T": ("Thr", "Treonina"),
    "W": ("Trp", "Triptófano"),
    "Y": ("Tyr", "Tirosina"),
    "V": ("Val", "Valina"),
    "*": ("Stop", "Parada"),
}

CODIGO_GENETICO = {
    "".join(codon): aa
    for codon, aa in zip(product(_BASES, repeat=3), _AMINOACIDOS_1)
}
CODON_INICIO = "AUG"
CODONES_PARADA = frozenset(c for c, aa in CODIGO_GENETICO.items() if aa == "*")


def info_codon(codon: str) -> dict:
    aa = CODIGO_GENETICO[codon]
    abrev, nombre = AMINOACIDOS[aa]
    return {
        "codon": codon,
        # Anticodón escrito 3'→5', alineado con el codón 5'→3'.
        "anticodon": complementar_arn(codon),
        "aa": aa,
        "abrev": abrev,
        "nombre": nombre,
        "parada": aa == "*",
    }


def traducir(arnm: str) -> dict:
    inicio = arnm.find(CODON_INICIO)
    base = {"arnm": arnm, "inicio": inicio, "codones": [], "proteina": "", "proteina_abrev": ""}

    if inicio == -1:
        return base | {
            "terminada": False,
            "fotos": [
                {
                    "etapa": "sin_inicio",
                    "titulo": "No hay codón de inicio",
                    "texto": (
                        "La subunidad menor del ribosoma recorre el ARNm 5'→3' y no encuentra "
                        "ningún AUG. Sin codón de inicio no empieza la traducción y no se "
                        "sintetiza ninguna proteína. Prueba con otra hebra molde o con "
                        "otra secuencia."
                    ),
                    "sitios": {},
                    "peptido": [],
                }
            ],
        }

    codones = []
    terminada = False
    i = inicio
    while i + 3 <= len(arnm):
        datos = info_codon(arnm[i : i + 3]) | {"posicion": i, "indice": len(codones)}
        codones.append(datos)
        i += 3
        if datos["parada"]:
            terminada = True
            break

    def sitios(e, p, a):
        valido = lambda x: x if x is not None and 0 <= x < len(codones) else None
        return {"E": valido(e), "P": valido(p), "A": valido(a)}

    peptido = [codones[0]["aa"]]
    fotos = [
        {
            "etapa": "iniciacion",
            "titulo": "Iniciación",
            "texto": (
                "La subunidad menor del ribosoma se une al extremo 5' del ARNm y lo recorre "
                f"hasta el primer codón AUG (nucleótido {inicio + 1}). Así queda fijado el "
                "marco de lectura. El ARNt iniciador, con anticodón 3'-UAC-5' y cargado con "
                "metionina, se sitúa en el sitio P. Después se acopla la subunidad mayor."
            ),
            "sitios": sitios(None, 0, 1),
            "peptido": list(peptido),
        }
    ]

    for j in range(1, len(codones)):
        c = codones[j]
        if c["parada"]:
            fotos.append(
                {
                    "etapa": "terminacion",
                    "titulo": f"Terminación en {c['codon']}",
                    "texto": (
                        f"El codón de parada {c['codon']} entra en el sitio A. Ningún ARNt lo "
                        "reconoce; en su lugar se une un factor de liberación. Este hace que "
                        "se hidrolice el enlace entre el último ARNt y la cadena, y la "
                        f"proteína de {len(peptido)} aminoácidos queda libre. Las subunidades "
                        "del ribosoma se separan."
                    ),
                    "sitios": sitios(None, j - 1, j),
                    "peptido": list(peptido),
                }
            )
            continue

        peptido.append(c["aa"])
        fotos.append(
            {
                "etapa": "elongacion",
                "titulo": f"Codón {j + 1}: {c['codon']} → {c['abrev']}",
                "texto": (
                    f"Un ARNt con anticodón 3'-{c['anticodon']}-5' se empareja con el codón "
                    f"{c['codon']} en el sitio A y aporta {c['nombre'].lower()}. La "
                    "peptidil transferasa (el ARNr de la subunidad mayor) forma el enlace "
                    "peptídico. Después el ribosoma avanza un codón (translocación): el ARNt "
                    "con la cadena pasa al sitio P y el ARNt vacío sale por el sitio E."
                ),
                "sitios": sitios(j - 1, j, j + 1),
                "peptido": list(peptido),
            }
        )

    if not terminada:
        fotos.append(
            {
                "etapa": "sin_parada",
                "titulo": "Fin del ARNm sin codón de parada",
                "texto": (
                    "El ribosoma llega al extremo 3' del ARNm sin encontrar un codón de "
                    "parada en este marco de lectura. La cadena de aminoácidos queda "
                    "incompleta. En la célula, estos ribosomas bloqueados se rescatan y "
                    "la proteína incompleta se degrada."
                ),
                "sitios": sitios(len(codones) - 2, len(codones) - 1, None),
                "peptido": list(peptido),
            }
        )

    proteina = "".join(peptido)
    return base | {
        "codones": codones,
        "terminada": terminada,
        "proteina": proteina,
        "proteina_abrev": "-".join(AMINOACIDOS[a][0] for a in proteina),
        "fotos": fotos,
    }
