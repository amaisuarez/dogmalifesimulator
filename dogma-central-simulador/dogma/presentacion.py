# -*- coding: utf-8 -*-
"""
dogma.presentacion
===================
Convierte los objetos de resultado del simulador en texto legible,
con formato de bloques por bases (para facilitar la lectura de
secuencias largas) y separadores claros entre etapas.
"""

from .simulador import ResultadoReplicacion, ResultadoTranscripcion, ResultadoTraduccion

ANCHO = 78


def _titulo(texto: str, relleno: str = "=") -> str:
    linea = relleno * ANCHO
    return f"\n{linea}\n{texto.center(ANCHO)}\n{linea}"


def _subtitulo(texto: str) -> str:
    return f"\n--- {texto} " + "-" * max(0, ANCHO - len(texto) - 5)


def formatear_en_bloques(seq: str, tam_bloque: int = 10) -> str:
    """Inserta un espacio cada `tam_bloque` bases, útil para leer ADN/ARN largos."""
    return " ".join(seq[i:i + tam_bloque] for i in range(0, len(seq), tam_bloque))


def mostrar_entrada(hebra_A: str, hebra_B: str) -> str:
    out = [_titulo("MOLÉCULA DE ADN DE PARTIDA (doble hélice)")]
    out.append(f"Hebra A 5'-{formatear_en_bloques(hebra_A)}-3'")
    out.append(f"Hebra B 3'-{formatear_en_bloques(hebra_B[::-1])}-5'   (complementaria y antiparalela)")
    return "\n".join(out)


def mostrar_replicacion(r: ResultadoReplicacion) -> str:
    out = [_titulo("1. REPLICACIÓN DEL ADN  (ADN -> ADN)")]
    for e in r.eventos:
        out.append(f"* {e}")

    out.append(_subtitulo("Hebra líder (síntesis continua)"))
    out.append(f"Molde   (Hebra A) 5'-{formatear_en_bloques(r.hebra_A_parental)}-3'")
    out.append(f"Nueva   (líder)   3'-{formatear_en_bloques(r.hebra_lider_nueva[::-1])}-5'")

    out.append(_subtitulo("Hebra retrasada (fragmentos de Okazaki)"))
    out.append(f"Molde (Hebra B, leída 3'->5' desde el extremo de la horquilla):")
    for f in r.fragmentos_okazaki:
        out.append(
            f"  Fragmento {f.indice:>2}: molde={f.posicion_molde:<10} "
            f"cebador={f.cebador_arn:<14} -> ADN nuevo={f.fragmento_adn}"
        )
    out.append(f"Hebra retrasada final (fragmentos ya unidos por la ligasa):")
    out.append(f"  5'-{formatear_en_bloques(r.hebra_retrasada_nueva)}-3'")

    out.append(_subtitulo("Moléculas de ADN hijas (semiconservativas)"))
    out.append("Dúplex 1: hebra parental A + hebra líder nueva")
    out.append(f"  5'-{formatear_en_bloques(r.duplex_1['hebra_1'])}-3'")
    out.append(f"  3'-{formatear_en_bloques(r.duplex_1['hebra_2'][::-1])}-5'")
    out.append("Dúplex 2: hebra parental B + hebra retrasada nueva")
    out.append(f"  3'-{formatear_en_bloques(r.duplex_2['hebra_1'][::-1])}-5'")
    out.append(f"  5'-{formatear_en_bloques(r.duplex_2['hebra_2'])}-3'")
    return "\n".join(out)


def mostrar_transcripcion(t: ResultadoTranscripcion) -> str:
    out = [_titulo("2. TRANSCRIPCIÓN  (ADN -> ARNm)")]
    for e in t.eventos:
        out.append(f"* {e}")
    out.append(_subtitulo("Complementariedad molde -> ARNm"))
    out.append(f"Hebra molde       3'-{formatear_en_bloques(t.hebra_molde[::-1])}-5'")
    out.append(f"ARNm sintetizado  5'-{formatear_en_bloques(t.arnm)}-3'")
    out.append(f"(Hebra codificante, para referencia: 5'-{formatear_en_bloques(t.hebra_codificante)}-3')")
    return "\n".join(out)


def mostrar_traduccion(tr: ResultadoTraduccion) -> str:
    out = [_titulo("3. TRADUCCIÓN  (ARNm -> proteína)")]
    for e in tr.eventos:
        out.append(f"* {e}")

    if not tr.proteina:
        out.append("\nNo se ha podido sintetizar ninguna proteína completa.")
        return "\n".join(out)

    out.append(_subtitulo("Lectura de codones"))
    codones_txt = " | ".join(tr.codones)
    out.append(f"Codones leídos: {codones_txt}")

    out.append(_subtitulo("Cadena polipeptídica resultante"))
    out.append("Secuencia (3 letras): " + "-".join(a.nombre3 for a in tr.proteina)
                + (f"-STOP({tr.codon_paro})" if tr.codon_paro else ""))
    out.append("Secuencia (1 letra):  " + "".join(a.letra1 for a in tr.proteina))
    out.append(f"Longitud: {len(tr.proteina)} aminoácidos")
    return "\n".join(out)


def mostrar_resumen_final(r: ResultadoReplicacion, t: ResultadoTranscripcion, tr: ResultadoTraduccion) -> str:
    out = [_titulo("RESUMEN DEL FLUJO DE INFORMACIÓN GENÉTICA", relleno="#")]
    out.append(f"ADN (hebra codificante) : 5'-{r.hebra_A_parental}-3'")
    out.append(f"        (replicación)   -> 2 moléculas de ADN hijas idénticas")
    out.append(f"ARNm  (transcripción)   : 5'-{t.arnm}-3'")
    if tr.proteina:
        out.append(f"Proteína (traducción)  : {tr.secuencia_1letra}  "
                    f"({'-'.join(a.nombre3 for a in tr.proteina)})")
    else:
        out.append("Proteína (traducción)  : no se generó (sin marco de lectura válido)")
    return "\n".join(out)
