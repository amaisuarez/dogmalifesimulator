"""Versión en texto del simulador.

Uso:
    python -m dogma ATGGCC...            # secuencia como argumento
    python -m dogma --aleatoria          # gen aleatorio
    python -m dogma SECUENCIA --molde superior --fragmento 6 --cebador 2
"""

import argparse
import sys

from .secuencias import generar_gen_aleatorio
from .simulador import simular

FLECHA = "\n" + " " * 4 + "↓\n"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m dogma", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("secuencia", nargs="?", help="hebra superior del ADN (5'→3')")
    ap.add_argument("--aleatoria", action="store_true", help="usar un gen aleatorio")
    ap.add_argument("--molde", choices=["inferior", "superior"], default="inferior")
    ap.add_argument("--fragmento", type=int, default=8, help="nt por fragmento de Okazaki")
    ap.add_argument("--cebador", type=int, default=3, help="nt por cebador de ARN")
    ap.add_argument("--pasos", action="store_true", help="mostrar cada paso de cada etapa")
    args = ap.parse_args(argv)

    if args.aleatoria or not args.secuencia:
        secuencia = generar_gen_aleatorio()
    else:
        secuencia = args.secuencia

    try:
        r = simular(secuencia, args.fragmento, args.cebador, args.molde)
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1

    rep, tr, td = r["replicacion"], r["transcripcion"], r["traduccion"]

    def pasos(fotos):
        if args.pasos:
            for i, f in enumerate(fotos, 1):
                print(f"  [{i}] {f['titulo']}\n      {f['texto']}")

    print("ADN PARENTAL")
    print(f"  5'-{rep['superior']}-3'")
    print(f"  3'-{rep['inferior']}-5'")
    print(FLECHA + "REPLICACIÓN (ADN → ADN)")
    pasos(rep["fotos"])
    for fr in rep["fragmentos"]:
        print(f"  Fragmento de Okazaki {fr['numero']}: posiciones {fr['inicio'] + 1}-{fr['fin']}"
              f" ({fr['longitud']} nt)")
    for h in rep["hijas"]:
        print(f"  {h['nombre']} (hebra parental: {h['parental']})")
        print(f"    5'-{h['superior']}-3'\n    3'-{h['inferior']}-5'")
    print(f"  Copia fiel: {'sí' if rep['copia_fiel'] else 'no'}")
    print(FLECHA + "TRANSCRIPCIÓN (ADN → ARN)")
    pasos(tr["fotos"])
    print(f"  Hebra molde: {tr['hebra_molde']}")
    print(f"  ARNm 5'-{tr['arnm']}-3'")
    print(FLECHA + "TRADUCCIÓN (ARN → proteína)")
    pasos(td["fotos"])
    if td["inicio"] == -1:
        print("  No se ha encontrado ningún codón AUG.")
        return 0
    print("  Codones: " + " ".join(c["codon"] for c in td["codones"]))
    print("  Anticod: " + " ".join(c["anticodon"] for c in td["codones"]))
    print(f"  Proteína: {td['proteina_abrev']}")
    print(f"            {td['proteina']} ({len(td['proteina'])} aa)")
    if not td["terminada"]:
        print("  Aviso: no se alcanzó ningún codón de parada.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
