"""Navegacion de la capa quimica en `nucleo.get_druggable_targets`.

`ddt` (Tanimoto) y `dds` (subestructura) unen CLUSTERS identitarios
(createDB/01-02). La semilla tiene que consultarlas con el `cluster_id` de la
droga, no con su `drug_id`. Hasta el 2026-09-24 se usaba el `drug_id`: como las
dos numeraciones se solapan no fallaba, devolvia los vecinos de otro cluster.

Tablas sinteticas: la droga 7 vive en el cluster 100 (Tanimoto) y 200
(subestructura); su vecino real es la droga 8 (cluster 101 / 201), con blanco
T_REAL. Existe ademas un cluster numerado 7 —el id de la droga— cuyo vecino es
la droga 9 con blanco T_SENUELO: si la semilla lo trae, se esta buscando con el
id equivocado.
"""
import sys, pathlib, unittest
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))   # comun/

import pandas as pd

import nucleo as nf


def tablas():
    tclus = pd.DataFrame({"drug": [7.0, 8.0, 9.0, 50.0], "cluster_id": [100, 101, 102, 7]})
    sclus = pd.DataFrame({"drug": [7, 8, 9, 50], "cluster_id": [200, 201, 202, 7]})
    ddt = pd.DataFrame({"clusID1": [100, 7], "clusID2": [101, 102], "weight": [0.9, 0.95]})
    dds = pd.DataFrame({"from": [200, 7], "to": [201, 202], "weight": [1.0, 1.0]}).set_index("from", drop=False)
    posdt = pd.DataFrame({"drug_id": [8, 8, 9, 9], "target_id": ["T_REAL", "T_REAL2", "T_SENUELO", "T_SENUELO2"]})
    st = pd.DataFrame({"target_id": ["T_REAL", "T_REAL2", "T_SENUELO", "T_SENUELO2"], "sp_id": ["x"] * 4})
    negdt = pd.DataFrame({"drug_id": [], "target_id": []})
    return dict(posdt=posdt, st=st, negdt=negdt, dds=dds, ddt=ddt, tclus=tclus, sclus=sclus)


class TestVecinosPorCluster(unittest.TestCase):

    def semilla(self):
        t = tablas()
        return nf.get_druggable_targets(doi=[7], use_neighbors=True, b_max_consolidate=False, **t)["7"]

    def test_trae_el_vecino_real(self):
        s = self.semilla()
        self.assertIsNotNone(s)
        self.assertIn("T_REAL", set(s.loc[s["type"] == "tani", "target_id"]))
        self.assertIn("T_REAL", set(s.loc[s["type"] == "sub", "target_id"]))

    def test_no_trae_el_cluster_con_el_id_de_la_droga(self):
        s = self.semilla()
        self.assertNotIn("T_SENUELO", set(s["target_id"]))

    def test_peso_es_la_similitud(self):
        s = self.semilla()
        w = s.loc[(s["type"] == "tani") & (s["target_id"] == "T_REAL"), "w"].iloc[0]
        self.assertAlmostEqual(w, 0.9)


if __name__ == "__main__":
    unittest.main()
