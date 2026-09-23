# -*- coding: utf-8 -*-
"""
dogma.biologia
==============
Tablas y constantes biológicas usadas por el simulador del dogma central:
- Reglas de complementariedad ADN-ADN y ADN-ARN.
- Código genético estándar (64 codones -> aminoácido).
- Codón de inicio y codones de parada.

Mantener estas tablas separadas del motor de simulación facilita su
verificación y su reutilización en los tres procesos (replicación,
transcripción y traducción).
"""

# ---------------------------------------------------------------------------
# Complementariedad de bases
# ---------------------------------------------------------------------------

# ADN molde -> ADN de nueva síntesis (replicación). A-T y G-C (o T-A y C-G).
COMPLEMENTO_ADN = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
}

# ADN molde -> ARN de nueva síntesis (transcripción). La Timina del molde
# se aparea con Adenina en el ARN, y la Adenina del molde con Uracilo
# (el ARN no contiene Timina).
COMPLEMENTO_ADN_A_ARN = {
    "A": "U",
    "T": "A",
    "G": "C",
    "C": "G",
}

BASES_ADN_VALIDAS = set(COMPLEMENTO_ADN.keys())

# ---------------------------------------------------------------------------
# Código genético estándar (64 codones de ARNm)
# Cada codón se mapea a: (nombre tres letras, letra IUPAC, tipo)
# tipo puede ser "aa" (aminoácido) o "stop" (codón de terminación)
# ---------------------------------------------------------------------------

CODIGO_GENETICO = {
    # Fenilalanina / Leucina
    "UUU": ("Phe", "F", "aa"), "UUC": ("Phe", "F", "aa"),
    "UUA": ("Leu", "L", "aa"), "UUG": ("Leu", "L", "aa"),
    "CUU": ("Leu", "L", "aa"), "CUC": ("Leu", "L", "aa"),
    "CUA": ("Leu", "L", "aa"), "CUG": ("Leu", "L", "aa"),
    # Isoleucina / Metionina (inicio)
    "AUU": ("Ile", "I", "aa"), "AUC": ("Ile", "I", "aa"), "AUA": ("Ile", "I", "aa"),
    "AUG": ("Met", "M", "start"),
    # Valina
    "GUU": ("Val", "V", "aa"), "GUC": ("Val", "V", "aa"),
    "GUA": ("Val", "V", "aa"), "GUG": ("Val", "V", "aa"),
    # Serina
    "UCU": ("Ser", "S", "aa"), "UCC": ("Ser", "S", "aa"),
    "UCA": ("Ser", "S", "aa"), "UCG": ("Ser", "S", "aa"),
    "AGU": ("Ser", "S", "aa"), "AGC": ("Ser", "S", "aa"),
    # Prolina
    "CCU": ("Pro", "P", "aa"), "CCC": ("Pro", "P", "aa"),
    "CCA": ("Pro", "P", "aa"), "CCG": ("Pro", "P", "aa"),
    # Treonina
    "ACU": ("Thr", "T", "aa"), "ACC": ("Thr", "T", "aa"),
    "ACA": ("Thr", "T", "aa"), "ACG": ("Thr", "T", "aa"),
    # Alanina
    "GCU": ("Ala", "A", "aa"), "GCC": ("Ala", "A", "aa"),
    "GCA": ("Ala", "A", "aa"), "GCG": ("Ala", "A", "aa"),
    # Tirosina / Stop
    "UAU": ("Tyr", "Y", "aa"), "UAC": ("Tyr", "Y", "aa"),
    "UAA": ("Ocre", "*", "stop"), "UAG": ("Ambar", "*", "stop"),
    # Histidina / Glutamina
    "CAU": ("His", "H", "aa"), "CAC": ("His", "H", "aa"),
    "CAA": ("Gln", "Q", "aa"), "CAG": ("Gln", "Q", "aa"),
    # Asparagina / Lisina
    "AAU": ("Asn", "N", "aa"), "AAC": ("Asn", "N", "aa"),
    "AAA": ("Lys", "K", "aa"), "AAG": ("Lys", "K", "aa"),
    # Ácido aspártico / Ácido glutámico
    "GAU": ("Asp", "D", "aa"), "GAC": ("Asp", "D", "aa"),
    "GAA": ("Glu", "E", "aa"), "GAG": ("Glu", "E", "aa"),
    # Cisteína / Stop / Triptófano
    "UGU": ("Cys", "C", "aa"), "UGC": ("Cys", "C", "aa"),
    "UGA": ("Opal", "*", "stop"),
    "UGG": ("Trp", "W", "aa"),
    # Arginina
    "CGU": ("Arg", "R", "aa"), "CGC": ("Arg", "R", "aa"),
    "CGA": ("Arg", "R", "aa"), "CGG": ("Arg", "R", "aa"),
    "AGA": ("Arg", "R", "aa"), "AGG": ("Arg", "R", "aa"),
    # Glicina
    "GGU": ("Gly", "G", "aa"), "GGC": ("Gly", "G", "aa"),
    "GGA": ("Gly", "G", "aa"), "GGG": ("Gly", "G", "aa"),
}

CODON_INICIO = "AUG"
CODONES_PARADA = {"UAA", "UAG", "UGA"}


def es_base_adn_valida(base: str) -> bool:
    return base.upper() in BASES_ADN_VALIDAS
