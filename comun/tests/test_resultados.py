"""Los resultados de la corrida, contra lo que se sabe de antemano.

A diferencia del resto de `comun/tests/`, estas pruebas necesitan que los
notebooks se hayan corrido: se saltean solas si no hay salidas. Son el paso de
`resultados/` de una carpeta de archivos a algo que afirma cosas verificables.
"""
import pathlib
import sys
import unittest

import pandas as pd

RAIZ = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(RAIZ / "comun"))
sys.path.insert(0, str(RAIZ / "genome_prioritization"))

import tdr as C  # noqa: E402

GP = C.SALIDAS / "genome_prioritization_out"
AN = C.SALIDAS / "analiceDB_out"
HAY_BARRIDO = (GP / "01_meta.json").exists()
HAY_DIM = (AN / "01_dimensiones_capas.csv").exists()


@unittest.skipUnless(HAY_BARRIDO, "sin barrido: correr genome_prioritization/01")
class TestControlExterno(unittest.TestCase):
    """El control independiente de `control/genoma_completo_v5/`."""

    @classmethod
    def setUpClass(cls):
        import funciones_genome as fg
        cls.fg = fg
        cls.optimos = fg.optimos_por_especie(fg.cargar_barrido(GP, "01"))
        cls.control = fg.comparar_con_control(cls.optimos)

    def test_estan_las_dieciseis_especies(self):
        self.assertEqual(len(self.optimos), 16, "el barrido no cubrio las 16 especies")

    def test_los_optimos_reproducen_el_control(self):
        """Los parametros optimos tienen que coincidir en las 16, y la AUC01
        exactamente: el control se calculo de forma independiente, sobre la base
        sin codificar. Si esto se rompe, el port dejo de ser fiel."""
        faltan = self.control.get("control", pd.Series(dtype=str)).eq("falta").sum()
        self.assertEqual(int(faltan), 0, "faltan especies en control/")
        self.assertEqual(int(self.control["coincide_params"].sum()), len(self.control),
                         "hay especies cuyo optimo no coincide con el control")
        self.assertLess(float(self.control["delta"].abs().max()), 1e-9,
                        "la AUC01 se movio respecto del control")

    def test_todas_superan_el_azar(self):
        self.assertTrue((self.optimos["AUC01"] > 0.5).all(),
                        "alguna especie no supera el azar en su optimo")


@unittest.skipUnless(HAY_DIM, "sin descriptivos: correr analiceDB/01")
class TestDimensiones(unittest.TestCase):
    """Lo que el recorte NO debia mover."""

    @classmethod
    def setUpClass(cls):
        d = pd.read_csv(AN / "01_dimensiones_capas.csv")
        cls.d = dict(zip(d["entidad"], d["n"]))

    def test_la_capa_de_anotaciones_no_se_movio(self):
        """El recorte saca compuestos, no proteinas ni anotaciones: estos tres
        numeros tienen que ser identicos a los de antes del recorte."""
        self.assertEqual(int(self.d["categorías InterPro"]), 14083)
        self.assertEqual(int(self.d["categorías OrthoMCL"]), 76157)
        self.assertEqual(int(self.d["enlaces proteína-categoría"]), 769661)

    def test_la_evidencia_positiva_sobrevivio_entera(self):
        """El criterio del recorte conserva toda componente con bioactividad
        positiva: los positivos y los blancos druggables no pueden bajar."""
        self.assertEqual(int(self.d["enlaces E_DP positivos"]), 248545)
        self.assertEqual(int(self.d["proteínas druggables"]), 19650)


class TestQuinasas(unittest.TestCase):
    """La marca de quinasa: es una revelacion parcial del diccionario privado
    (PROVENANCE.md §3.15), asi que conviene que no se mueva sin que se note."""

    def test_la_constante_no_se_movio(self):
        self.assertEqual(len(C.ANN_QUINASA), 391)
        self.assertEqual(len(C.anotaciones_quinasa()), 391)

    def test_todo_accession_es_de_interpro(self):
        """El valor de cada entrada es un accession `IPR` de siete digitos: si
        alguien pega un nombre de dominio en vez del accession, se ve aca."""
        for acc in C.ANN_QUINASA.values():
            self.assertRegex(acc, r"^IPR\d{6,7}$")

    def test_los_codigos_estan_en_la_capa_de_anotaciones(self):
        """Las 391 tienen que existir en `DB/` y ser de tipo Domain: si no, la
        marca no alcanzaria ninguna proteina y el analisis de beta daria vacio
        sin lanzar nada."""
        ip = pd.read_csv(C.RAIZ / "DB" / "04a_interpro.csv")
        ip = ip[ip["type"] == C.IPR_DOMAIN]
        sel = ip[ip["ann"].isin(set(C.ANN_QUINASA))]
        self.assertEqual(sel["ann"].nunique(), 391)
        self.assertEqual(sel["target_id"].nunique(), 13921)


@unittest.skipUnless((GP / "03_meta.json").exists(),
                     "sin comparacion de beta: correr genome_prioritization/03")
class TestComparacionBeta(unittest.TestCase):
    """El contraste de `genome_prioritization/03` tiene que ser a parametros
    controlados: lo unico que cambia entre las dos corridas es beta."""

    @classmethod
    def setUpClass(cls):
        cls.params = pd.read_csv(GP / "03_params_comparacion.csv")
        cls.resumen = pd.read_csv(GP / "03_resumen_beta.csv")

    def test_solo_se_movio_beta(self):
        self.assertEqual(len(self.params), 16)
        self.assertTrue((self.params["beta"] == 1.0).all())

    def test_la_auc01_de_beta_alto_sale_del_barrido(self):
        """La AUC01 que el notebook recalcula con beta = 1 tiene que ser la
        misma que el barrido del notebook 01 ya habia medido para esa
        combinacion: si difieren, el contexto o la metrica cambiaron."""
        j = self.params.merge(self.resumen, on="especie")
        self.assertEqual(len(j), 16)
        for _, r in j.iterrows():
            self.assertAlmostEqual(r["AUC01"], r["auc01_b1"], places=9,
                                   msg=f"especie {r['especie']}")

    def test_el_conteo_del_top_tolera_empates(self):
        """Con beta = 0 los scores se degeneran y el corte del top-100 cae
        dentro de un bloque empatado; el conteo es el esperado bajo desempate
        uniforme, asi que nunca puede pasarse de las 100 posiciones ni del
        total de quinasas de la especie."""
        for col in ("quinasas_top100_b1", "quinasas_top100_b0"):
            self.assertTrue((self.resumen[col] <= 100).all())
            self.assertTrue((self.resumen[col] <= self.resumen["n_quinasas"]).all())
            self.assertTrue((self.resumen[col] >= 0).all())


if __name__ == "__main__":
    unittest.main()
