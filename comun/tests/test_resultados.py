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


if __name__ == "__main__":
    unittest.main()
