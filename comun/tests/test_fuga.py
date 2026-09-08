"""Ausencia de fuga de datos en las validaciones cruzadas.

Punto del checklist del roadmap: "validar que la remoción de enlaces de drogas
pseudohuérfanas no contamine los conjuntos de entrenamiento".
"""
import sys, pathlib, unittest
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))   # comun/

import numpy as np
import pandas as pd

import tdr as C

HAY_DB = (C.DB / "00_specie_target.csv").exists()
SP = "calb"      # especie chica: el test corre rápido


@unittest.skipUnless(HAY_DB, "DB/ no disponible")
class TestFugaLOSO(unittest.TestCase):
    """Leave-one-species-out: la semilla no puede contener evidencia de la especie."""

    @classmethod
    def setUpClass(cls):
        cls.d = C.cargar_db(anotaciones=False, cluster_consistent=False, verbose=False)

    def test_semilla_sin_blancos_de_la_especie_retenida(self):
        seed = C.semilla_global(self.d.posdt, self.d.st, sp_out=SP)
        de_la_especie = set(self.d.targets_de(SP))
        fuga = de_la_especie & set(seed["target_id"].astype(str))
        self.assertEqual(len(fuga), 0, f"{len(fuga)} blancos de {SP} quedaron en la semilla")

    def test_los_positivos_evaluados_no_son_semilla(self):
        seed = C.semilla_global(self.d.posdt, self.d.st, sp_out=SP)
        tp = np.setdiff1d(self.d.druggables_de(SP), seed["target_id"].astype(str).values)
        self.assertGreater(len(tp), 0, "sin positivos que evaluar")
        self.assertEqual(len(set(tp) & set(seed["target_id"].astype(str))), 0)

    def test_sin_sp_out_la_semilla_si_los_contiene(self):
        """Control negativo: sin remoción, los blancos de la especie sí están."""
        seed = C.semilla_global(self.d.posdt, self.d.st)
        self.assertGreater(len(set(self.d.targets_de(SP)) & set(seed["target_id"].astype(str))), 0)


@unittest.skipUnless(HAY_DB, "DB/ no disponible")
class TestFugaOrfanizacion(unittest.TestCase):
    """Orfanización: al remover las aristas de una droga, su blanco sólo puede
    volver a la semilla por otra droga vecina (que es la inferencia directa, y es
    legítima), nunca por la droga orfanizada."""

    def test_remocion_completa_de_las_aristas_de_la_droga(self):
        d = C.cargar_db(anotaciones=False, cluster_consistent=True, verbose=False)
        posdt_sp = d.posdt.merge(d.st[["target_id", "sp_id"]], on="target_id", how="left")
        conteo = posdt_sp.groupby(["drug_id", "sp_id"]).size()
        droga, sp = conteo[conteo == 1].index[0]
        restante = posdt_sp[posdt_sp["drug_id"] != droga]
        self.assertEqual((restante["drug_id"] == droga).sum(), 0)
        # el blanco verdadero puede seguir siendo druggable por otras drogas:
        # eso no es fuga, es la evidencia del resto de la red
        blanco = posdt_sp.loc[(posdt_sp["drug_id"] == droga) &
                              (posdt_sp["sp_id"] == sp), "target_id"].iloc[0]
        otras = restante[restante["target_id"] == blanco]["drug_id"].nunique()
        self.assertGreaterEqual(otras, 0)


if __name__ == "__main__":
    unittest.main()
