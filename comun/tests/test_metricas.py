"""Tests de las métricas: AUC parcial y corrección de McClish.

Punto del checklist del roadmap: "Correr tests unitarios sobre el cálculo de la
AUC01 corregida por McClish" — un clasificador al azar debe dar 0.5 y el
clasificador perfecto 1.0.
"""
import sys, pathlib, unittest
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))   # comun/

import numpy as np

import tdr as C


class TestMcClish(unittest.TestCase):

    def test_azar_da_medio(self):
        """Un clasificador sin información: la pAUC normalizada tiende a FPRmax/2."""
        rng = np.random.default_rng(0)
        n = 40_000
        y = (rng.random(n) < 0.05).astype(int)
        s = rng.random(n)                       # score independiente de la etiqueta
        self.assertAlmostEqual(C.auc01_mcclish(y, s), 0.5, delta=0.03)

    def test_perfecto_da_uno(self):
        y = np.r_[np.ones(200), np.zeros(2000)].astype(int)
        s = np.r_[np.ones(200), np.zeros(2000)] + np.linspace(0, 0.1, 2200)
        self.assertAlmostEqual(C.auc01_mcclish(y, s), 1.0, places=6)

    def test_escala_conocida(self):
        """La transformación mapea el mínimo a 0.5 y el máximo a 1."""
        self.assertAlmostEqual(C.mcclish(0.05), 0.5)
        self.assertAlmostEqual(C.mcclish(1.0), 1.0)
        self.assertAlmostEqual(C.mcclish(0.525), 0.75)

    def test_monotona(self):
        """La corrección no cambia el orden: los óptimos del barrido se preservan."""
        p = np.linspace(0.05, 1.0, 50)
        m = C.mcclish(p)
        self.assertTrue(np.all(np.diff(m) > 0))

    def test_formula_invertida_del_notebook(self):
        """La fórmula de la versión v3 (0.5*(1-(p-0.005)/0.095)) es decreciente y
        cae fuera de [0.5, 1]: este test documenta por qué se reemplaza."""
        p = np.array([0.05, 0.5, 1.0])
        mala = 0.5 * (1 - (p - 0.005) / 0.095)
        self.assertTrue(np.all(np.diff(mala) < 0))          # decreciente
        self.assertTrue(mala.min() < 0)                      # sale del rango
        buena = C.mcclish(p)
        self.assertTrue(np.all((buena >= 0.5 - 1e-9) & (buena <= 1 + 1e-9)))

    def test_pauc_respeta_el_corte(self):
        """La pAUC sólo mira FPR <= 0.1: un cambio en la cola no la altera."""
        rng = np.random.default_rng(1)
        n = 5000
        y = (rng.random(n) < 0.1).astype(int)
        s = rng.random(n) + y * 0.5
        a = C.pauc_normalizada(y, s, 0.10)
        s2 = s.copy()
        peores = np.argsort(s2)[:1000]                       # reordena sólo la cola baja
        s2[peores] = s2[peores][::-1]
        b = C.pauc_normalizada(y, s2, 0.10)
        self.assertAlmostEqual(a, b, places=6)

    def test_bootstrap(self):
        rng = np.random.default_rng(2)
        n = 2000
        y = (rng.random(n) < 0.1).astype(int)
        s = rng.random(n) + y * 0.4
        b = C.bootstrap_auc01(y, s, n_boot=50)
        self.assertEqual(len(b), 50)
        self.assertTrue(np.nanmin(b) >= 0.4)


if __name__ == "__main__":
    unittest.main()
