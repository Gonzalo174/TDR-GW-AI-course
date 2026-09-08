"""Equivalencia entre la propagación de `nds_fun.nds()` y la proyección del paper.

El roadmap especifica la red proyectada como

    M^PP = M^bip S (M^bip)^T ,      S_ll = RS_l * k_l^(-beta)
    PS_i = sum_j w_j e_ji

mientras que `nds()` la implementa como una propagación en dos pasos con un
parámetro extra `lambda`. Este test verifica numéricamente la correspondencia
sobre una red de juguete:

    nds(lambda=0)_i * k_i  ==  (M S M^T w)_i

es decir, con lambda = 0 la propagación reproduce exactamente la proyección del
paper salvo el factor de normalización por grado de la proteína receptora, que
`nds()` aplica como k_i^(lambda-1).
"""
import sys, pathlib, unittest
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))   # comun/

import numpy as np
import pandas as pd

import tdr as C


def red_juguete():
    """6 proteínas, 4 categorías, pesos de semilla arbitrarios."""
    enlaces = [("P1", "IPR1"), ("P1", "OG1"), ("P2", "IPR1"), ("P3", "IPR1"),
               ("P3", "OG2"), ("P4", "OG1"), ("P4", "IPR2"), ("P5", "IPR2"),
               ("P6", "OG2"), ("P6", "IPR1")]
    sta = pd.DataFrame(enlaces, columns=["target_id", "ann"])
    sta["db"] = np.where(sta["ann"].str.startswith("IPR"), "ip", "omcl")
    sta["sp_id"] = "toy"
    seed = pd.DataFrame({"target_id": ["P1", "P2", "P4"], "w": [1.0, 2.0, 0.5]})
    cat_rs = pd.DataFrame({"ann": ["IPR1", "IPR2", "OG1", "OG2"],
                           "RS": [1.0, 0.4, 0.7, 0.2]})
    return sta, seed, cat_rs


def proyeccion_explicita(sta, seed, cat_rs, beta):
    prot = sorted(sta["target_id"].unique())
    cats = sorted(sta["ann"].unique())
    M = np.zeros((len(prot), len(cats)))
    for _, r in sta.iterrows():
        M[prot.index(r["target_id"]), cats.index(r["ann"])] = 1.0
    k_cat = M.sum(axis=0)
    rs = cat_rs.set_index("ann").loc[cats, "RS"].values
    S = np.diag(rs * k_cat ** (-beta))
    MPP = M @ S @ M.T
    w = np.zeros(len(prot))
    for _, r in seed.iterrows():
        w[prot.index(r["target_id"])] = r["w"]
    return prot, MPP @ w


class TestProyeccion(unittest.TestCase):

    def test_equivalencia_con_lambda_cero(self):
        sta, seed, cat_rs = red_juguete()
        for beta in (0.0, 1.0, 1.5):
            with self.subTest(beta=beta):
                prot, ps = proyeccion_explicita(sta, seed, cat_rs, beta)
                r = C.propagar(sta, seed, cat_rs, beta=beta, lambda_=0.0)
                r = r.set_index("target_id").loc[prot]
                k_prot = sta.groupby("target_id").size().loc[prot].values
                np.testing.assert_allclose(r["score"].values * k_prot, ps, rtol=1e-10)

    def test_beta_penaliza_categorias_grandes(self):
        """beta = 1 (G'rk) debe bajar el peso relativo de la categoría más grande."""
        sta, seed, cat_rs = red_juguete()
        r0 = C.propagar(sta, seed, cat_rs, beta=0.0, lambda_=1.0).set_index("target_id")
        r1 = C.propagar(sta, seed, cat_rs, beta=1.0, lambda_=1.0).set_index("target_id")
        # P2 sólo cuelga de IPR1, la categoría más grande (4 proteínas)
        self.assertLess(r1.loc["P2", "score"] / r1["score"].sum(),
                        r0.loc["P2", "score"] / r0["score"].sum())

    def test_todas_las_proteinas_reciben_score(self):
        sta, seed, cat_rs = red_juguete()
        r = C.propagar(sta, seed, cat_rs, beta=1.0, lambda_=1.0)
        self.assertEqual(len(r), sta["target_id"].nunique())
        self.assertTrue((r["score"] >= 0).all())


if __name__ == "__main__":
    unittest.main()
