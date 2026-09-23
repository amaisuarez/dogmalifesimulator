# -*- coding: utf-8 -*-
"""
dogma.simulador
================
Motor de simulación del dogma central de la biología molecular:

    ADN --(replicación)--> ADN
    ADN --(transcripción)--> ARNm
    ARNm --(traducción)--> Proteína

El objetivo es didáctico: cada método no solo calcula el resultado, sino
que registra ("narra") los pasos intermedios -moléculas, enzimas y reglas
de complementariedad empleadas- en una lista de eventos que la capa de
presentación (main.py) convierte en texto legible por el usuario.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Optional

from .biologia import (
    COMPLEMENTO_ADN,
    COMPLEMENTO_ADN_A_ARN,
    CODIGO_GENETICO,
    CODON_INICIO,
    CODONES_PARADA,
    es_base_adn_valida,
)


# ---------------------------------------------------------------------------
# Utilidades de secuencias
# ---------------------------------------------------------------------------

def limpiar_secuencia(seq: str) -> str:
    """Elimina espacios/saltos de línea y pasa a mayúsculas."""
    return "".join(seq.split()).upper()


def complemento_adn(seq: str) -> str:
    return "".join(COMPLEMENTO_ADN[b] for b in seq)


def reverso_complemento_adn(seq: str) -> str:
    return complemento_adn(seq)[::-1]


def validar_adn(seq: str) -> None:
    invalidas = sorted({b for b in seq if not es_base_adn_valida(b)})
    if invalidas:
        raise ValueError(
            "La secuencia de ADN contiene símbolos no válidos: "
            f"{', '.join(invalidas)}. Solo se admiten A, T, G, C."
        )
    if len(seq) == 0:
        raise ValueError("La secuencia de ADN no puede estar vacía.")


# ---------------------------------------------------------------------------
# Estructuras de resultado
# ---------------------------------------------------------------------------

@dataclass
class FragmentoOkazaki:
    indice: int
    cebador_arn: str
    fragmento_adn: str
    posicion_molde: str  # tramo del molde (hebra retrasada) copiado


@dataclass
class ResultadoReplicacion:
    hebra_A_parental: str
    hebra_B_parental: str
    cebador_lider: str
    hebra_lider_nueva: str
    fragmentos_okazaki: List[FragmentoOkazaki]
    hebra_retrasada_nueva: str
    duplex_1: Dict[str, str]
    duplex_2: Dict[str, str]
    eventos: List[str] = field(default_factory=list)


@dataclass
class ResultadoTranscripcion:
    hebra_molde: str
    hebra_codificante: str
    arnm: str
    eventos: List[str] = field(default_factory=list)


@dataclass
class Aminoacido:
    codon: str
    nombre3: str
    letra1: str


@dataclass
class ResultadoTraduccion:
    arnm_leido: str
    marco_lectura_inicio: int
    codones: List[str]
    proteina: List[Aminoacido]
    codon_paro: Optional[str]
    eventos: List[str] = field(default_factory=list)

    @property
    def secuencia_3letras(self) -> str:
        return "-".join(a.nombre3 for a in self.proteina)

    @property
    def secuencia_1letra(self) -> str:
        return "".join(a.letra1 for a in self.proteina)


# ---------------------------------------------------------------------------
# Simulador principal
# ---------------------------------------------------------------------------

class SimuladorDogmaCentral:
    """
    Orquesta los tres procesos del dogma central sobre una molécula de
    ADN de doble hebra proporcionada por el usuario (hebra codificante,
    5'->3').
    """

    def __init__(self, hebra_codificante_5_3: str):
        seq = limpiar_secuencia(hebra_codificante_5_3)
        validar_adn(seq)
        self.hebra_A = seq                          # 5' -> 3' (dada por el usuario)
        self.hebra_B = reverso_complemento_adn(seq)  # 5' -> 3' (complementaria, antiparalela)

    # ------------------------------------------------------------------
    # 1) REPLICACIÓN DEL ADN
    # ------------------------------------------------------------------
    def replicar(self, tamano_fragmento_okazaki: int = 8) -> ResultadoReplicacion:
        """
        Simula la replicación semiconservativa de la doble hélice.

        Simplificación pedagógica (se indica explícitamente): se asume que
        la horquilla de replicación avanza en un único sentido a lo largo
        del fragmento de ADN proporcionado.
          - La Hebra A actúa como molde de la HEBRA LÍDER: se sintetiza de
            forma continua, en dirección 5'->3', tras la apertura de la
            doble hélice por la helicasa y la colocación de un único
            cebador de ARN por la primasa.
          - La Hebra B actúa como molde de la HEBRA RETRASADA: al ser
            copiada en sentido contrario al avance de la horquilla, se
            sintetiza de forma discontinua en fragmentos de Okazoki, cada
            uno iniciado por su propio cebador de ARN (primasa) y
            polimerizado 5'->3' (ADN polimerasa III). La ADN ligasa une
            después los fragmentos, tras la eliminación de los cebadores.
        """
        eventos: List[str] = []
        eventos.append(
            "La helicasa reconoce el origen de replicación y rompe los "
            "puentes de hidrógeno entre las bases, abriendo la doble "
            "hélice y generando una horquilla de replicación."
        )
        eventos.append(
            "Las proteínas de unión a ADN monocatenario (SSB) estabilizan "
            "las hebras separadas para impedir que vuelvan a aparearse, "
            "mientras la topoisomerasa alivia la tensión torsional por "
            "delante de la horquilla."
        )

        # --- Hebra líder (síntesis continua, molde = Hebra A) ---
        cebador_lider = "PPP-" + complemento_adn(self.hebra_A[:5]).replace("T", "U") + "..."
        eventos.append(
            "PRIMASA: sintetiza un único cebador de ARN complementario al "
            f"inicio de la Hebra A (molde), permitiendo a la ADN "
            "polimerasa III comenzar la síntesis."
        )
        # La hebra líder recién sintetizada es, base a base, la
        # complementaria de la Hebra A (molde).
        hebra_lider_nueva = complemento_adn(self.hebra_A)
        eventos.append(
            "ADN POLIMERASA III sintetiza la HEBRA LÍDER de forma "
            "continua en dirección 5'->3', leyendo la Hebra A molde en "
            "dirección 3'->5' y aplicando la complementariedad A-T / G-C."
        )

        # --- Hebra retrasada (síntesis discontinua, molde = Hebra B) ---
        eventos.append(
            "Como la Hebra B se lee en el mismo sentido que avanza la "
            "horquilla, no puede copiarse de forma continua: la ADN "
            "polimerasa III sintetiza la HEBRA RETRASADA en fragmentos "
            "cortos y discontinuos (FRAGMENTOS DE OKAZAKI), cada uno "
            "iniciado por un cebador de ARN independiente."
        )
        molde_retrasada_invertido = self.hebra_B[::-1]  # se procesa desde el extremo 3' de B
        fragmentos: List[FragmentoOkazaki] = []
        n = tamano_fragmento_okazaki
        trozos = [
            molde_retrasada_invertido[i:i + n]
            for i in range(0, len(molde_retrasada_invertido), n)
        ]
        hebra_retrasada_partes = []
        for idx, trozo in enumerate(trozos, start=1):
            cebador = complemento_adn(trozo[: min(3, len(trozo))]).replace("T", "U")
            fragmento_adn = complemento_adn(trozo)
            hebra_retrasada_partes.append(fragmento_adn)
            fragmentos.append(
                FragmentoOkazaki(
                    indice=idx,
                    cebador_arn=f"5'-{cebador}-3' (ARN)",
                    fragmento_adn=fragmento_adn,
                    posicion_molde=trozo,
                )
            )
            eventos.append(
                f"  Fragmento de Okazaki #{idx}: primasa coloca cebador de "
                f"ARN -> ADN polimerasa III extiende {len(fragmento_adn)} "
                "nt de ADN sobre el molde (Hebra B)."
            )

        hebra_retrasada_nueva = "".join(hebra_retrasada_partes)
        eventos.append(
            "RNasa H / ADN polimerasa I retiran los cebadores de ARN y "
            "los sustituyen por ADN; a continuación, la ADN LIGASA sella "
            "los huecos entre fragmentos de Okazaki consecutivos, "
            "formando una hebra retrasada continua."
        )

        duplex_1 = {"hebra_1": self.hebra_A, "hebra_2": hebra_lider_nueva}
        duplex_2 = {"hebra_1": self.hebra_B, "hebra_2": hebra_retrasada_nueva}
        eventos.append(
            "Resultado: dos moléculas de ADN hijas, cada una formada por "
            "una hebra parental original y una hebra recién sintetizada "
            "(replicación SEMICONSERVATIVA), ambas idénticas a la "
            "molécula de ADN original."
        )

        return ResultadoReplicacion(
            hebra_A_parental=self.hebra_A,
            hebra_B_parental=self.hebra_B,
            cebador_lider=cebador_lider,
            hebra_lider_nueva=hebra_lider_nueva,
            fragmentos_okazaki=fragmentos,
            hebra_retrasada_nueva=hebra_retrasada_nueva,
            duplex_1=duplex_1,
            duplex_2=duplex_2,
            eventos=eventos,
        )

    # ------------------------------------------------------------------
    # 2) TRANSCRIPCIÓN (ADN -> ARNm)
    # ------------------------------------------------------------------
    def transcribir(self, usar_hebra: str = "B") -> ResultadoTranscripcion:
        """
        Sintetiza el ARN mensajero a partir de la hebra molde indicada.

        usar_hebra="B" (por defecto): la Hebra B actúa como MOLDE y la
        Hebra A es entonces la hebra CODIFICANTE (su secuencia coincide
        con el ARNm, salvo T->U). Esta es la convención más habitual en
        los libros de texto introductorios.
        """
        eventos: List[str] = []
        if usar_hebra.upper() == "B":
            molde = self.hebra_B
            codificante = self.hebra_A
        elif usar_hebra.upper() == "A":
            molde = self.hebra_A
            codificante = self.hebra_B
        else:
            raise ValueError("usar_hebra debe ser 'A' o 'B'")

        eventos.append(
            "La ARN POLIMERASA reconoce una región promotora en el ADN y "
            "se une a la doble hélice, provocando su apertura local."
        )
        eventos.append(
            f"Se selecciona como HEBRA MOLDE la hebra {'B' if usar_hebra.upper()=='B' else 'A'} "
            "(se lee en dirección 3'->5'); la hebra complementaria actúa "
            "como HEBRA CODIFICANTE (su secuencia es idéntica al ARNm "
            "resultante, sustituyendo T por U)."
        )

        # El ARNm se sintetiza 5'->3' leyendo el molde 3'->5'. Como el
        # molde se ha almacenado en sentido 5'->3', se recorre invertido.
        arnm = "".join(COMPLEMENTO_ADN_A_ARN[b] for b in molde[::-1])
        eventos.append(
            "La ARN polimerasa sintetiza el ARNm en dirección 5'->3', "
            "aplicando complementariedad de bases sobre el molde "
            "(A(molde)->U, T(molde)->A, G(molde)->C, C(molde)->G)."
        )
        eventos.append(
            "Al alcanzar una secuencia terminadora, la ARN polimerasa se "
            "libera del ADN junto con el transcrito de ARNm recién "
            "sintetizado."
        )

        return ResultadoTranscripcion(
            hebra_molde=molde,
            hebra_codificante=codificante,
            arnm=arnm,
            eventos=eventos,
        )

    # ------------------------------------------------------------------
    # 3) TRADUCCIÓN (ARNm -> proteína)
    # ------------------------------------------------------------------
    def traducir(self, arnm: str) -> ResultadoTraduccion:
        """
        Traduce el ARNm en una secuencia de aminoácidos, leyendo codones
        a partir del primer codón de inicio (AUG) encontrado, hasta un
        codón de parada.
        """
        eventos: List[str] = []
        arnm = arnm.upper().replace(" ", "")

        pos_inicio = arnm.find(CODON_INICIO)
        if pos_inicio == -1:
            eventos.append(
                "No se ha encontrado ningún codón de inicio (AUG) en el "
                "ARNm: la traducción no puede comenzar."
            )
            return ResultadoTraduccion(
                arnm_leido=arnm,
                marco_lectura_inicio=-1,
                codones=[],
                proteina=[],
                codon_paro=None,
                eventos=eventos,
            )

        eventos.append(
            f"El RIBOSOMA (subunidad pequeña) se une al ARNm y localiza el "
            f"primer codón de inicio AUG en la posición {pos_inicio + 1}, "
            "que fija el marco de lectura."
        )
        eventos.append(
            "Un ARN de transferencia iniciador (ARNt-Met) reconoce el "
            "codón AUG por complementariedad con su anticodón; se ensambla "
            "la subunidad ribosómica grande, formando el ribosoma completo."
        )

        codones: List[str] = []
        proteina: List[Aminoacido] = []
        codon_paro: Optional[str] = None

        i = pos_inicio
        while i + 3 <= len(arnm):
            codon = arnm[i:i + 3]
            codones.append(codon)
            info = CODIGO_GENETICO.get(codon)
            if info is None:
                eventos.append(f"Codón {codon} no reconocido; se detiene la traducción.")
                break
            nombre3, letra1, tipo = info
            if tipo == "stop":
                codon_paro = codon
                eventos.append(
                    f"Codón de PARADA {codon} reconocido por un factor de "
                    "liberación (no existe ARNt complementario): el "
                    "ribosoma libera la cadena polipeptídica y se "
                    "disocia del ARNm."
                )
                break
            proteina.append(Aminoacido(codon=codon, nombre3=nombre3, letra1=letra1))
            eventos.append(
                f"Codón {codon} -> ARNt con anticodón complementario aporta "
                f"{nombre3} ({letra1}); la peptidil-transferasa del "
                "ribosoma forma el enlace peptídico con la cadena en "
                "crecimiento."
            )
            i += 3

        if codon_paro is None and proteina:
            eventos.append(
                "Se ha alcanzado el final del ARNm sin encontrar un codón "
                "de parada explícito (secuencia incompleta)."
            )

        return ResultadoTraduccion(
            arnm_leido=arnm,
            marco_lectura_inicio=pos_inicio,
            codones=codones,
            proteina=proteina,
            codon_paro=codon_paro,
            eventos=eventos,
        )

    # ------------------------------------------------------------------
    # Flujo completo
    # ------------------------------------------------------------------
    def ejecutar_dogma_completo(
        self, tamano_fragmento_okazaki: int = 8, usar_hebra_molde: str = "B"
    ):
        replicacion = self.replicar(tamano_fragmento_okazaki)
        transcripcion = self.transcribir(usar_hebra_molde)
        traduccion = self.traducir(transcripcion.arnm)
        return replicacion, transcripcion, traduccion
