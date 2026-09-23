# -*- coding: utf-8 -*-
"""
Pruebas unitarias básicas del simulador del dogma central.
Ejecutar con:  python -m pytest tests/  (o  python -m unittest tests/test_simulador.py)
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from dogma.simulador import SimuladorDogmaCentral, reverso_complemento_adn, validar_adn


class TestComplementariedad(unittest.TestCase):
    def test_reverso_complemento(self):
        self.assertEqual(reverso_complemento_adn("ATGC"), "GCAT")
        self.assertEqual(reverso_complemento_adn("AAAA"), "TTTT")

    def test_validar_adn_rechaza_caracteres_invalidos(self):
        with self.assertRaises(ValueError):
            validar_adn("ATGX")

    def test_validar_adn_rechaza_vacia(self):
        with self.assertRaises(ValueError):
            validar_adn("")


class TestReplicacion(unittest.TestCase):
    def setUp(self):
        self.sim = SimuladorDogmaCentral("ATGGCTAAGCGTTTCTGGTAA")

    def test_hebra_B_es_complementaria_antiparalela(self):
        self.assertEqual(self.sim.hebra_B, reverso_complemento_adn(self.sim.hebra_A))

    def test_duplex_hijos_reconstruyen_original(self):
        r = self.sim.replicar(tamano_fragmento_okazaki=5)
        # Dúplex 1: hebra A parental + hebra líder nueva (complementaria de A)
        self.assertEqual(r.duplex_1["hebra_1"], self.sim.hebra_A)
        self.assertEqual(r.duplex_1["hebra_2"], r.hebra_lider_nueva)
        # La hebra retrasada, tras unir los fragmentos, debe reconstruir la Hebra A
        self.assertEqual(r.hebra_retrasada_nueva, self.sim.hebra_A)

    def test_fragmentos_okazaki_generados(self):
        r = self.sim.replicar(tamano_fragmento_okazaki=5)
        self.assertGreater(len(r.fragmentos_okazaki), 1)
        total_bases = sum(len(f.fragmento_adn) for f in r.fragmentos_okazaki)
        self.assertEqual(total_bases, len(self.sim.hebra_B))


class TestTranscripcion(unittest.TestCase):
    def test_arnm_no_contiene_timina(self):
        sim = SimuladorDogmaCentral("ATGGCTAAGCGTTTCTGGTAA")
        t = sim.transcribir()
        self.assertNotIn("T", t.arnm)

    def test_longitud_arnm_igual_a_adn(self):
        sim = SimuladorDogmaCentral("ATGGCTAAGCGTTTCTGGTAA")
        t = sim.transcribir()
        self.assertEqual(len(t.arnm), len(sim.hebra_A))


class TestTraduccion(unittest.TestCase):
    def test_traduccion_correcta_con_stop(self):
        sim = SimuladorDogmaCentral("ATGGCTAAGCGTTTCTGGTAA")
        t = sim.transcribir()
        tr = sim.traducir(t.arnm)
        self.assertEqual(tr.secuencia_1letra, "MAKRFW")
        self.assertEqual(tr.codon_paro, "UAA")

    def test_sin_codon_inicio(self):
        sim = SimuladorDogmaCentral("GGGCCCTTTAAA")
        t = sim.transcribir()
        tr = sim.traducir(t.arnm)
        self.assertEqual(len(tr.proteina), 0)
        self.assertEqual(tr.marco_lectura_inicio, -1)


if __name__ == "__main__":
    unittest.main()
