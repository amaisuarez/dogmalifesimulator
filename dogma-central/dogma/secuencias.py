"""Utilidades básicas para trabajar con secuencias de ácidos nucleicos.

Aquí se definen las reglas de complementariedad de bases que usan
el resto de módulos:

* ADN - ADN (replicación):   A-T, T-A, C-G, G-C
* ADN - ARN (transcripción): A-U, T-A, C-G, G-C
* ARN - ARN (codón-anticodón en la traducción): A-U, U-A, C-G, G-C
"""

from __future__ import annotations

import random

BASES_ADN = frozenset("ATCG")
BASES_ARN = frozenset("AUCG")

# Emparejamiento de Watson y Crick entre dos hebras de ADN.
COMPLEMENTO_ADN = {"A": "T", "T": "A", "C": "G", "G": "C"}

# Base de ARN que la ARN polimerasa (o la primasa) coloca frente a
# cada base de la hebra molde de ADN. El uracilo sustituye a la timina.
MOLDE_ADN_A_ARN = {"A": "U", "T": "A", "C": "G", "G": "C"}

# Emparejamiento entre bases de ARN (codón del ARNm con anticodón del ARNt).
COMPLEMENTO_ARN = {"A": "U", "U": "A", "C": "G", "G": "C"}

CODONES_PARADA_ADN = ("TAA", "TAG", "TGA")


def limpiar_secuencia(texto: str) -> str:
    """Normaliza una secuencia de ADN introducida por el usuario.

    Acepta formato FASTA (ignora las líneas que empiezan por '>'),
    mayúsculas o minúsculas, espacios, saltos de línea y números de
    posición. Lanza ``ValueError`` si hay caracteres que no son bases de ADN.
    """
    if texto is None:
        raise ValueError("No se ha introducido ninguna secuencia.")

    lineas = [l for l in texto.splitlines() if not l.strip().startswith(">")]
    limpia = "".join(c for c in "".join(lineas).upper() if not c.isspace() and not c.isdigit())

    if not limpia:
        raise ValueError("La secuencia está vacía.")

    invalidos = sorted(set(limpia) - BASES_ADN)
    if invalidos:
        if "U" in invalidos:
            raise ValueError(
                "La secuencia contiene uracilo (U). El simulador parte de una "
                "molécula de ADN: usa T en lugar de U."
            )
        raise ValueError(
            "Caracteres no válidos en la secuencia: " + ", ".join(invalidos)
            + ". Solo se admiten A, T, C y G."
        )
    return limpia


def complementar(secuencia: str) -> str:
    """Hebra complementaria de ADN, alineada posición a posición (sin invertir)."""
    return "".join(COMPLEMENTO_ADN[b] for b in secuencia)


def complementaria_inversa(secuencia: str) -> str:
    """Hebra complementaria leída en sentido 5'→3'."""
    return complementar(secuencia)[::-1]


def arn_desde_molde(molde: str) -> str:
    """Bases de ARN complementarias a una hebra molde de ADN (alineadas)."""
    return "".join(MOLDE_ADN_A_ARN[b] for b in molde)


def complementar_arn(secuencia: str) -> str:
    """Complementaria de una secuencia de ARN (alineada, sin invertir)."""
    return "".join(COMPLEMENTO_ARN[b] for b in secuencia)


def contenido_gc(secuencia: str) -> float:
    """Porcentaje de G + C en la secuencia."""
    if not secuencia:
        return 0.0
    gc = sum(1 for b in secuencia if b in "GC")
    return round(100 * gc / len(secuencia), 1)


def generar_gen_aleatorio(codones: int = 12, flanco: int = 8, semilla: int | None = None) -> str:
    """Genera una secuencia de ADN que contiene un marco de lectura abierto.

    La secuencia tiene la forma: flanco 5' + ATG + codones + codón de
    parada + flanco 3'. Los codones internos nunca son de parada, y los
    flancos se eligen de forma que no introduzcan un ATG antes del real.
    """
    rng = random.Random(semilla)
    codones = max(1, min(codones, 80))
    flanco = max(0, min(flanco, 30))

    def aleatoria(n: int) -> str:
        return "".join(rng.choice("ATCG") for _ in range(n))

    # Flanco 5' sin ATG (para que el primer AUG del ARNm sea el del gen).
    cinco = aleatoria(flanco)
    while "ATG" in cinco:
        cinco = aleatoria(flanco)

    cuerpo = []
    while len(cuerpo) < codones:
        codon = aleatoria(3)
        if codon not in CODONES_PARADA_ADN:
            cuerpo.append(codon)

    parada = rng.choice(CODONES_PARADA_ADN)
    tres = aleatoria(flanco)
    return cinco + "ATG" + "".join(cuerpo) + parada + tres
