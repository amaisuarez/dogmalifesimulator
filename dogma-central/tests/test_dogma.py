import pytest

from dogma import replicar, simular, traducir, transcribir
from dogma.secuencias import (
    complementar,
    complementaria_inversa,
    generar_gen_aleatorio,
    limpiar_secuencia,
)
from dogma.traduccion import CODIGO_GENETICO

HBB = "CTGACGCCACCATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCTAAGCTT"


# ------------------------------------------------------------- secuencias
def test_complementariedad_adn():
    assert complementar("ATCG") == "TAGC"
    assert complementaria_inversa("AACG") == "CGTT"


def test_limpiar_acepta_fasta_minusculas_y_espacios():
    assert limpiar_secuencia(">gen\natg gcc\n12 taa") == "ATGGCCTAA"


@pytest.mark.parametrize("entrada", ["", "ATGX", "AUGC"])
def test_limpiar_rechaza_entradas_invalidas(entrada):
    with pytest.raises(ValueError):
        limpiar_secuencia(entrada)


def test_gen_aleatorio_tiene_marco_completo():
    for semilla in range(30):
        adn = generar_gen_aleatorio(codones=10, semilla=semilla)
        prot = traducir(adn.replace("T", "U"))
        assert prot["terminada"]
        assert len(prot["proteina"]) == 11  # Met + 10 codones


# ------------------------------------------------------------- replicación
@pytest.mark.parametrize("frag,ceb", [(8, 3), (5, 2), (4, 1), (20, 10)])
def test_replicacion_produce_copias_fieles(frag, ceb):
    r = replicar(HBB, tam_fragmento=frag, tam_cebador=ceb)
    assert r["copia_fiel"]
    for hija in r["hijas"]:
        assert complementar(hija["superior"]) == hija["inferior"]


def test_replicacion_es_semiconservativa():
    r = replicar(HBB)
    h1, h2 = r["hijas"]
    assert h1["superior"] == HBB and h1["parental"] == "superior"
    assert h2["inferior"] == complementar(HBB) and h2["parental"] == "inferior"


def test_fragmentos_de_okazaki_cubren_la_rezagada():
    r = replicar(HBB, tam_fragmento=8, tam_cebador=3)
    frags = r["fragmentos"]
    assert frags[0]["inicio"] == 0 and frags[-1]["fin"] == len(HBB)
    for a, b in zip(frags, frags[1:]):
        assert a["fin"] == b["inicio"]
    assert len(frags) == -(-len(HBB) // 8)


def test_cebadores_son_de_arn_y_se_eliminan():
    r = replicar(HBB, tam_fragmento=8, tam_cebador=3)
    fotos = r["fotos"]
    ceb = next(f for f in fotos if f["etapa"] == "cebador_lider")
    assert ceb["lider"]["tipos"][:3] == "rrr"
    assert set(ceb["lider"]["bases"][:3]) <= set("AUCG")
    assert "U" in "".join(f["rezagada"]["bases"] for f in fotos) or "U" in ceb["lider"]["bases"]
    final = fotos[-1]
    assert set(final["lider"]["tipos"]) == {"d"}
    assert set(final["rezagada"]["tipos"]) == {"d"}
    assert final["mellas"] == []


def test_cadena_lider_usa_un_solo_cebador():
    r = replicar(HBB, tam_fragmento=6, tam_cebador=2)
    assert sum(1 for f in r["fotos"] if f["etapa"] == "cebador_lider") == 1
    assert sum(1 for f in r["fotos"] if f["etapa"] == "cebador_rezagada") == len(r["fragmentos"])


def test_intervienen_todas_las_enzimas():
    r = replicar(HBB)
    claves = {e["clave"] for f in r["fotos"] for e in f["enzimas"]}
    assert {"helicasa", "topoisomerasa", "ssb", "primasa", "pol3", "pol1", "ligasa"} <= claves


def test_parametros_invalidos():
    with pytest.raises(ValueError):
        replicar(HBB, tam_fragmento=4, tam_cebador=4)


# ------------------------------------------------------------ transcripción
def test_transcripcion_con_molde_inferior():
    t = transcribir("ATGCCT", "inferior")
    assert t["arnm"] == "AUGCCU"
    assert t["molde"] == "TACGGA"


def test_transcripcion_con_molde_superior():
    t = transcribir("ATGCCT", "superior")
    assert t["arnm"] == "AGGCAU"  # complementaria inversa de ATGCCT, con U


def test_arnm_es_complementario_al_molde():
    t = transcribir(HBB, "inferior")
    pares = {"A": "U", "T": "A", "C": "G", "G": "C"}
    assert "".join(pares[b] for b in t["molde"]) == t["arnm"]
    assert "T" not in t["arnm"]


# --------------------------------------------------------------- traducción
def test_codigo_genetico_estandar():
    assert len(CODIGO_GENETICO) == 64
    assert CODIGO_GENETICO["AUG"] == "M"
    assert CODIGO_GENETICO["UGG"] == "W"
    assert {c for c, a in CODIGO_GENETICO.items() if a == "*"} == {"UAA", "UAG", "UGA"}


def test_traduccion_hasta_codon_de_parada():
    t = traducir("GGAUGUUUUGGUAAGCC")
    assert t["inicio"] == 2
    assert t["proteina"] == "MFW"
    assert t["terminada"]
    assert [c["codon"] for c in t["codones"]] == ["AUG", "UUU", "UGG", "UAA"]


def test_anticodones():
    t = traducir("AUGGCC")
    assert [c["anticodon"] for c in t["codones"]] == ["UAC", "CGG"]


def test_sin_codon_de_inicio():
    t = traducir("CCCGGGUUU")
    assert t["inicio"] == -1 and t["proteina"] == ""


def test_sin_codon_de_parada():
    t = traducir("AUGGCCGCC")
    assert t["proteina"] == "MAA" and not t["terminada"]


# ----------------------------------------------------------- flujo completo
def test_flujo_completo_beta_globina():
    r = simular(HBB)
    assert r["traduccion"]["proteina"] == "MVHLTPEEKSA"


def test_mutacion_falciforme_cambia_glu_por_val():
    mutada = HBB.replace("CCTGAGGAG", "CCTGTGGAG")
    assert simular(mutada)["traduccion"]["proteina"] == "MVHLTPVEKSA"


def test_gen_en_la_hebra_inferior():
    r = simular("GGATCCTTATTTGCTGGACCAGTAGAACTTCATGGTGGCGAATTC", hebra_molde="superior")
    assert r["traduccion"]["proteina"] == "MKFYWSSK"


def test_limites_de_longitud():
    with pytest.raises(ValueError):
        simular("ATG")
    with pytest.raises(ValueError):
        simular("A" * 301)


def test_api_flask():
    from app import app

    cliente = app.test_client()
    ok = cliente.post("/api/simular", json={"secuencia": HBB})
    assert ok.status_code == 200 and ok.get_json()["traduccion"]["proteina"] == "MVHLTPEEKSA"
    mal = cliente.post("/api/simular", json={"secuencia": "HOLA"})
    assert mal.status_code == 400 and "error" in mal.get_json()
    assert cliente.get("/").status_code == 200
