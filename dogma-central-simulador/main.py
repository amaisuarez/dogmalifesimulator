#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simulador del Dogma Central de la Biología Molecular
=====================================================
Práctica 1 - Bioinformática

Uso:
    python main.py                                  # modo interactivo, secuencia de ejemplo
    python main.py --secuencia ATGGCC...TAA          # simular una secuencia propia
    python main.py --archivo mi_adn.txt              # leer la secuencia de un fichero
    python main.py --okazaki 6 --sin-diagrama        # opciones adicionales

Ejecuta, en orden, los tres procesos del dogma central sobre una molécula
de ADN de doble hebra:
    1) Replicación del ADN (hebra líder + hebra retrasada/fragmentos de Okazaki)
    2) Transcripción (ADN -> ARNm)
    3) Traducción (ARNm -> proteína)

y muestra el resultado paso a paso por consola, además de generar (por
defecto) una imagen PNG con el resumen visual del proceso completo.
"""

import argparse
import sys
import os

from dogma.simulador import SimuladorDogmaCentral
from dogma import presentacion
from dogma.diagrama import generar_diagrama

SECUENCIA_EJEMPLO = "CGTACGTACGATGGCTAAGCGTTTCTGGTAATTTTAAACCC"


def leer_secuencia_interactiva() -> str:
    print(presentacion._titulo("SIMULADOR DEL DOGMA CENTRAL — MODO INTERACTIVO", "="))
    print(
        "\nIntroduce una secuencia de ADN (hebra 5'->3', solo bases A, T, G, C)."
        "\nPulsa ENTER sin escribir nada para usar una secuencia de ejemplo.\n"
    )
    try:
        entrada = input("Secuencia de ADN > ").strip()
    except EOFError:
        entrada = ""
    return entrada if entrada else SECUENCIA_EJEMPLO


def construir_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Simulador del dogma central de la biología molecular "
                    "(replicación, transcripción y traducción).",
    )
    grupo_secuencia = parser.add_mutually_exclusive_group()
    grupo_secuencia.add_argument(
        "--secuencia", "-s", type=str, default=None,
        help="Secuencia de ADN (hebra codificante 5'->3'), p. ej. ATGGCCATT...",
    )
    grupo_secuencia.add_argument(
        "--archivo", "-f", type=str, default=None,
        help="Ruta a un fichero de texto que contiene la secuencia de ADN.",
    )
    parser.add_argument(
        "--okazaki", "-k", type=int, default=8,
        help="Tamaño (en nucleótidos) de cada fragmento de Okazaki simulado (por defecto: 8).",
    )
    parser.add_argument(
        "--hebra-molde", choices=["A", "B"], default="B",
        help="Hebra usada como molde en la transcripción (por defecto: B).",
    )
    parser.add_argument(
        "--sin-diagrama", action="store_true",
        help="No generar la imagen PNG resumen del proceso.",
    )
    parser.add_argument(
        "--salida", "-o", type=str, default="resultados",
        help="Carpeta donde guardar el diagrama y el informe de texto (por defecto: ./resultados).",
    )
    return parser


def main() -> int:
    parser = construir_parser()
    args = parser.parse_args()

    # --- Obtener la secuencia de ADN ---
    if args.secuencia:
        secuencia = args.secuencia
    elif args.archivo:
        try:
            with open(args.archivo, "r", encoding="utf-8") as f:
                secuencia = f.read()
        except OSError as exc:
            print(f"No se pudo leer el archivo '{args.archivo}': {exc}", file=sys.stderr)
            return 1
    else:
        secuencia = leer_secuencia_interactiva()

    # --- Ejecutar el simulador ---
    try:
        simulador = SimuladorDogmaCentral(secuencia)
    except ValueError as exc:
        print(f"\nError en la secuencia de ADN: {exc}", file=sys.stderr)
        return 1

    replicacion, transcripcion, traduccion = simulador.ejecutar_dogma_completo(
        tamano_fragmento_okazaki=args.okazaki,
        usar_hebra_molde=args.hebra_molde,
    )

    # --- Construir el informe de texto ---
    bloques = [
        presentacion.mostrar_entrada(simulador.hebra_A, simulador.hebra_B),
        presentacion.mostrar_replicacion(replicacion),
        presentacion.mostrar_transcripcion(transcripcion),
        presentacion.mostrar_traduccion(traduccion),
        presentacion.mostrar_resumen_final(replicacion, transcripcion, traduccion),
    ]
    informe = "\n".join(bloques) + "\n"
    print(informe)

    # --- Guardar salidas ---
    os.makedirs(args.salida, exist_ok=True)
    ruta_informe = os.path.join(args.salida, "informe_simulacion.txt")
    with open(ruta_informe, "w", encoding="utf-8") as f:
        f.write(informe)
    print(f"\n[Informe de texto guardado en: {ruta_informe}]")

    if not args.sin_diagrama:
        ruta_png = os.path.join(args.salida, "diagrama_dogma_central.png")
        generar_diagrama(replicacion, transcripcion, traduccion, ruta_png)
        print(f"[Diagrama guardado en: {ruta_png}]")

    return 0


if __name__ == "__main__":
    sys.exit(main())
