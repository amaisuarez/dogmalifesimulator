"""Simulación de la replicación semiconservativa del ADN (ADN → ADN).

Modelo empleado
---------------
La molécula se dibuja con la hebra superior en sentido 5'→3' (de izquierda
a derecha) y la inferior en sentido 3'→5'. Hay una única horquilla que parte
de un origen situado en el extremo izquierdo y avanza hacia la derecha.

* Las ADN polimerasas solo sintetizan en sentido 5'→3', añadiendo
  nucleótidos al extremo 3'-OH de una cadena ya existente.
* La hebra inferior (3'→5' en el sentido de avance) sirve de molde a la
  **cadena líder**, que crece de forma continua hacia la horquilla y solo
  necesita un cebador.
* La hebra superior (5'→3' en el sentido de avance) sirve de molde a la
  **cadena rezagada**, que tiene que crecer alejándose de la horquilla.
  Por eso se sintetiza a trozos, los **fragmentos de Okazaki**, cada uno con
  su propio cebador de ARN.

Cada paso se guarda como una "foto" del estado, que la interfaz reproduce.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

from .secuencias import COMPLEMENTO_ADN, MOLDE_ADN_A_ARN, complementar

VACIO, CEBADOR, ADN = "-", "r", "d"


@dataclass
class Fragmento:
    """Un fragmento de Okazaki de la cadena rezagada (posiciones [inicio, fin))."""

    numero: int
    inicio: int
    fin: int
    cebador_inicio: int
    cebador_fin: int

    @property
    def longitud(self) -> int:
        return self.fin - self.inicio


def _enzima(clave: str, inicio: int, fin: int, hebra: str) -> dict:
    return {"clave": clave, "inicio": inicio, "fin": fin, "hebra": hebra}


class _Replicador:
    def __init__(self, superior: str, tam_fragmento: int, tam_cebador: int):
        self.sup = superior
        self.inf = complementar(superior)
        self.n = len(superior)
        self.f = tam_fragmento
        self.p = tam_cebador

        # Hebras nuevas: base, tipo (vacío / cebador de ARN / ADN) y fragmento.
        self.lider_b = ["·"] * self.n
        self.lider_t = [VACIO] * self.n
        self.rez_b = ["·"] * self.n
        self.rez_t = [VACIO] * self.n
        self.rez_frag = [0] * self.n

        self.horquilla = 0
        self.mellas: set[int] = set()  # fronteras i-1 | i sin enlace fosfodiéster
        self.fragmentos: list[Fragmento] = []
        self.fotos: list[dict] = []

    # ------------------------------------------------------------------ fotos
    def _enzimas_horquilla(self) -> list[dict]:
        h, n = self.horquilla, self.n
        lista = []
        if h < n:
            lista.append(_enzima("helicasa", h, h + 1, "horquilla"))
            if h + 2 < n:
                lista.append(_enzima("topoisomerasa", h + 2, min(n, h + 4), "horquilla"))
        return lista

    def _ssb(self, hebra: str, inicio: int, fin: int) -> list[dict]:
        return [_enzima("ssb", inicio, fin, hebra)] if fin > inicio else []

    def foto(self, etapa: str, titulo: str, texto: str, enzimas: list[dict]) -> None:
        self.fotos.append(
            {
                "etapa": etapa,
                "titulo": titulo,
                "texto": texto,
                "horquilla": self.horquilla,
                "lider": {"bases": "".join(self.lider_b), "tipos": "".join(self.lider_t)},
                "rezagada": {
                    "bases": "".join(self.rez_b),
                    "tipos": "".join(self.rez_t),
                    "fragmentos": list(self.rez_frag),
                },
                "mellas": sorted(self.mellas),
                "enzimas": self._enzimas_horquilla() + enzimas,
            }
        )

    # ------------------------------------------------------------ operaciones
    def _cebador_lider(self, inicio: int, fin: int) -> None:
        for i in range(inicio, fin):
            self.lider_b[i] = MOLDE_ADN_A_ARN[self.inf[i]]
            self.lider_t[i] = CEBADOR

    def _adn_lider(self, inicio: int, fin: int) -> None:
        for i in range(inicio, fin):
            self.lider_b[i] = COMPLEMENTO_ADN[self.inf[i]]
            self.lider_t[i] = ADN

    def _cebador_rezagada(self, inicio: int, fin: int, frag: int) -> None:
        for i in range(inicio, fin):
            self.rez_b[i] = MOLDE_ADN_A_ARN[self.sup[i]]
            self.rez_t[i] = CEBADOR
            self.rez_frag[i] = frag

    def _adn_rezagada(self, inicio: int, fin: int, frag: int) -> None:
        for i in range(inicio, fin):
            self.rez_b[i] = COMPLEMENTO_ADN[self.sup[i]]
            self.rez_t[i] = ADN
            self.rez_frag[i] = frag

    # --------------------------------------------------------------- proceso
    def ejecutar(self) -> None:
        n, f, p = self.n, self.f, self.p

        self.foto(
            "origen",
            "Reconocimiento del origen de replicación",
            "Las proteínas iniciadoras reconocen el origen de replicación (aquí, el "
            "extremo izquierdo de la molécula) y reclutan a la helicasa. La doble "
            "hélice parental todavía está cerrada: cada base está unida a su "
            "complementaria por puentes de hidrógeno (A=T con 2, C≡G con 3).",
            [_enzima("iniciadoras", 0, min(n, 3), "horquilla")],
        )

        previo = 0
        lider_fin = 0
        numero = 0
        while previo < n:
            nuevo = min(n, previo + f)
            self.horquilla = nuevo
            numero += 1

            # 1. Apertura de la doble hélice.
            self.foto(
                "apertura",
                f"Apertura de la doble hélice (posiciones {previo + 1}–{nuevo})",
                "La helicasa avanza rompiendo los puentes de hidrógeno y separa las dos "
                "hebras parentales. Por delante, la topoisomerasa corta y vuelve a unir "
                "el ADN para aliviar el superenrollamiento. Las proteínas SSB se unen al "
                "ADN de cadena sencilla para que no vuelva a aparearse.",
                self._ssb("lider", lider_fin, nuevo) + self._ssb("rezagada", previo, nuevo),
            )

            # 2. Cebador único de la cadena líder.
            if numero == 1:
                fin_ceb = min(p, n)
                self._cebador_lider(0, fin_ceb)
                lider_fin = fin_ceb
                ceb = "".join(self.lider_b[:fin_ceb])
                self.foto(
                    "cebador_lider",
                    "Cebador de la cadena líder",
                    f"La primasa (una ARN polimerasa) sintetiza un cebador de ARN, 5'-{ceb}-3', "
                    "complementario a la hebra inferior. La ADN polimerasa no puede empezar "
                    "una cadena desde cero: necesita el extremo 3'-OH que le aporta el "
                    "cebador. La cadena líder solo necesita este cebador.",
                    [_enzima("primasa", 0, fin_ceb, "lider")]
                    + self._ssb("lider", fin_ceb, nuevo)
                    + self._ssb("rezagada", previo, nuevo),
                )

            # 3. Síntesis continua de la cadena líder.
            if lider_fin < nuevo:
                inicio = lider_fin
                self._adn_lider(inicio, nuevo)
                lider_fin = nuevo
                nuevas = "".join(self.lider_b[inicio:nuevo])
                self.foto(
                    "lider",
                    "Síntesis continua de la cadena líder",
                    f"La ADN polimerasa III lee la hebra inferior en sentido 3'→5' y añade "
                    f"desoxirribonucleótidos al extremo 3' de la cadena líder (5'…{nuevas}-3'). "
                    "Crece en el mismo sentido en que avanza la horquilla, así que la "
                    "síntesis es continua.",
                    [_enzima("pol3", nuevo - 1, nuevo, "lider")]
                    + self._ssb("rezagada", previo, nuevo),
                )

            # 4. Cebador del nuevo fragmento de Okazaki.
            largo = nuevo - previo
            ceb_ini = nuevo - min(p, largo)
            self._cebador_rezagada(ceb_ini, nuevo, numero)
            frag = Fragmento(numero, previo, nuevo, ceb_ini, nuevo)
            ceb = "".join(self.rez_b[ceb_ini:nuevo])[::-1]
            self.foto(
                "cebador_rezagada",
                f"Cebador del fragmento de Okazaki {numero}",
                "La hebra superior se lee 3'→5' en sentido contrario al avance de la "
                "horquilla, así que la nueva cadena tiene que crecer alejándose de ella. "
                f"La primasa coloca un cebador de ARN junto a la horquilla, 5'-{ceb}-3' "
                "(leído de derecha a izquierda).",
                [_enzima("primasa", ceb_ini, nuevo, "rezagada")]
                + self._ssb("rezagada", previo, ceb_ini),
            )

            # 5. Extensión del fragmento de Okazaki.
            if ceb_ini > previo:
                self._adn_rezagada(previo, ceb_ini, numero)
                destino = (
                    "el origen de replicación"
                    if numero == 1
                    else f"el cebador del fragmento {numero - 1}"
                )
                self.foto(
                    "okazaki",
                    f"Síntesis del fragmento de Okazaki {numero}",
                    "La ADN polimerasa III extiende el cebador en sentido 5'→3' (de derecha a "
                    f"izquierda en el dibujo) hasta llegar a {destino}. Resultado: un "
                    f"fragmento de {largo} nt ({largo - (nuevo - ceb_ini)} de ADN y "
                    f"{nuevo - ceb_ini} de ARN).",
                    [_enzima("pol3", previo, previo + 1, "rezagada")],
                )

            # 6-7. Retirada del cebador anterior y unión de fragmentos.
            if self.fragmentos:
                ant = self.fragmentos[-1]
                self._adn_rezagada(ant.cebador_inicio, ant.cebador_fin, numero)
                self.mellas.add(ant.cebador_inicio)
                self.foto(
                    "pol1",
                    f"Retirada del cebador del fragmento {ant.numero}",
                    "La ADN polimerasa I usa su actividad exonucleasa 5'→3' para retirar el "
                    f"cebador de ARN del fragmento {ant.numero}. A la vez, rellena el hueco "
                    f"con ADN extendiendo el extremo 3' del fragmento {numero}. Queda una "
                    "mella: dos nucleótidos contiguos sin enlace fosfodiéster.",
                    [_enzima("pol1", ant.cebador_inicio, ant.cebador_fin, "rezagada")],
                )
                self.mellas.discard(ant.cebador_inicio)
                inicio_lig = max(0, ant.cebador_inicio - 1)
                self.foto(
                    "ligasa",
                    f"Unión de los fragmentos {ant.numero} y {numero}",
                    "La ADN ligasa sella la mella formando el enlace fosfodiéster que "
                    "faltaba (consume ATP o NAD⁺). Los dos fragmentos de Okazaki quedan "
                    "unidos en una cadena continua.",
                    [_enzima("ligasa", inicio_lig, min(n, inicio_lig + 2), "rezagada")],
                )

            self.fragmentos.append(frag)
            previo = nuevo

        # 8. Cebadores de los extremos.
        fin_ceb_lider = min(p, n)
        ultimo = self.fragmentos[-1]
        self._adn_lider(0, fin_ceb_lider)
        self._adn_rezagada(ultimo.cebador_inicio, ultimo.cebador_fin, ultimo.numero)
        self.foto(
            "extremos",
            "Sustitución de los cebadores de los extremos",
            "Quedan dos cebadores de ARN: el de la cadena líder, en el origen, y el del "
            "último fragmento de Okazaki. En una molécula circular, o cuando hay "
            "replicación bidireccional, la horquilla vecina aporta el extremo 3'-OH que "
            "la ADN polimerasa I necesita para sustituirlos por ADN, y la ligasa sella "
            "las mellas. En el extremo de un cromosoma lineal, el hueco del último "
            "cebador no se puede rellenar; es el problema de la replicación de los "
            "extremos, que la telomerasa compensa en los telómeros. El simulador rellena "
            "ambos huecos para mostrar las moléculas hijas completas.",
            [
                _enzima("pol1", 0, fin_ceb_lider, "lider"),
                _enzima("pol1", ultimo.cebador_inicio, ultimo.cebador_fin, "rezagada"),
            ],
        )
        self.foto(
            "resultado",
            "Dos moléculas hijas idénticas",
            "La replicación es semiconservativa: cada molécula hija tiene una hebra "
            "parental y una hebra nueva. Además, la ADN polimerasa corrige errores "
            "con su actividad exonucleasa 3'→5' (corrección de pruebas), por lo que "
            "las dos copias son idénticas a la molécula original.",
            [],
        )


def replicar(superior: str, tam_fragmento: int = 8, tam_cebador: int = 3) -> dict:
    """Simula la replicación de la molécula cuya hebra superior (5'→3') se indica.

    ``tam_fragmento`` es la longitud de cada fragmento de Okazaki y
    ``tam_cebador`` la de cada cebador de ARN, ambos en nucleótidos. En la
    célula los fragmentos miden ~1000-2000 nt en bacterias y ~100-200 nt en
    eucariotas, y los cebadores ~10 nt. Aquí se usan valores pequeños para
    que el proceso quepa en pantalla.
    """
    if tam_fragmento < 2:
        raise ValueError("Los fragmentos de Okazaki deben medir al menos 2 nt.")
    if not 1 <= tam_cebador < tam_fragmento:
        raise ValueError("El cebador debe medir al menos 1 nt y ser más corto que el fragmento.")

    r = _Replicador(superior, tam_fragmento, tam_cebador)
    r.ejecutar()

    nueva_rezagada = "".join(r.rez_b)
    nueva_lider = "".join(r.lider_b)
    hijas = [
        {
            "nombre": "Molécula hija 1",
            "superior": superior,
            "inferior": nueva_rezagada,
            "parental": "superior",
        },
        {
            "nombre": "Molécula hija 2",
            "superior": nueva_lider,
            "inferior": r.inf,
            "parental": "inferior",
        },
    ]
    return {
        "superior": superior,
        "inferior": r.inf,
        "tam_fragmento": tam_fragmento,
        "tam_cebador": tam_cebador,
        "fragmentos": [asdict(fr) | {"longitud": fr.longitud} for fr in r.fragmentos],
        "hijas": hijas,
        "copia_fiel": nueva_lider == superior and nueva_rezagada == r.inf,
        "fotos": r.fotos,
    }
