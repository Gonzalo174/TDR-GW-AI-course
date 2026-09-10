"""Las figuras se dibujan desde las tablas de `resultados/`, nunca desde la base.

Es la regla que permite regenerar todas las figuras en segundos
(`python comun/figuras.py`) sin volver a pagar la corrida. Se rompe sin que nadie
lo note en cuanto alguien vuelve a dibujar dentro de una celda de corrida, o una
figura empieza a leer algo que no declaró: estas pruebas lo atrapan.

Las que dibujan se saltean solas para las figuras cuyas tablas todavía no están.
"""
import json
import pathlib
import re
import sys
import unittest
from unittest import mock

import pandas as pd

RAIZ = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(RAIZ / "comun"))

import figuras  # noqa: E402
import tdr  # noqa: E402

REG = figuras.cargar_registro()
NOTEBOOKS = sorted(RAIZ.glob("*/[0-9][0-9]_*.ipynb"))


def _celdas_de_codigo(nb):
    j = json.loads(nb.read_text())
    return ["".join(c["source"]) for c in j["cells"] if c["cell_type"] == "code"]


class TestRegistro(unittest.TestCase):

    def test_hay_figuras_de_los_cuatro_analisis(self):
        self.assertEqual({f.analisis for f in REG.values()}, set(figuras.MODULOS))

    def test_los_nombres_llevan_el_numero_del_notebook(self):
        for nombre in REG:
            with self.subTest(figura=nombre):
                self.assertRegex(nombre, r"^\d\d_f\d\d_[a-z0-9_]+$")

    def test_cada_figura_declara_lo_que_lee(self):
        for nombre, f in REG.items():
            with self.subTest(figura=nombre):
                self.assertTrue(f.lee, "sin tablas declaradas")
                for t in f.lee:
                    self.assertTrue(t.endswith(".csv"), t)

    def test_cada_figura_dice_que_la_verifica(self):
        for nombre, f in REG.items():
            with self.subTest(figura=nombre):
                self.assertTrue(f.verifica.strip(), "falta `verifica=` en @figura")

    def test_el_indice_esta_al_dia(self):
        """FIGURAS.md se genera del registro: si alguien agrega o cambia una
        figura sin regenerarlo, el índice que lee el corrector miente."""
        actual = (RAIZ / "FIGURAS.md").read_text()
        self.assertEqual(actual, figuras.indice(),
                         "FIGURAS.md desactualizado: correr `python comun/figuras.py --indice`")

    def test_los_notebooks_no_dibujan_por_su_cuenta(self):
        """Ninguna celda de notebook arma o guarda una figura: eso vive en
        `figuras_*.py`. Si vuelve a pasar, esa figura sólo se regenera corriendo
        el notebook entero."""
        prohibido = re.compile(r"plt\.subplots|tdr\.guardar\(|\.savefig\(|\.fig_[a-z]")
        for nb in NOTEBOOKS:
            for i, s in enumerate(_celdas_de_codigo(nb)):
                with self.subTest(notebook=nb.name, celda=i):
                    self.assertIsNone(prohibido.search(s), f"{nb.name}, celda {i}")

    def test_lo_que_dibujan_los_notebooks_esta_registrado(self):
        for nb in NOTEBOOKS:
            for s in _celdas_de_codigo(nb):
                for nombre in re.findall(r'F\.dibujar\("([^"]+)"\)', s):
                    with self.subTest(notebook=nb.name, figura=nombre):
                        self.assertIn(nombre, REG)

    def test_las_figuras_de_los_entregables_estan_registradas(self):
        """El informe y la página sólo pueden usar figuras que se regeneran
        desde las tablas."""
        tex = "".join((RAIZ / "informe" / f).read_text()
                      for f in ("informe.tex", "presentacion.tex"))
        pagina = (RAIZ / "informe" / "generar_pagina.py").read_text()
        usadas = set(re.findall(r"figuras/([A-Za-z0-9_]+)\.pdf", tex))
        usadas |= set(re.findall(r'"(\d\d_f\d\d_[a-z0-9_]+)\.png"', pagina))
        self.assertTrue(usadas)
        for nombre in sorted(usadas):
            with self.subTest(figura=nombre):
                self.assertIn(nombre, REG)


class TestDibujar(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.plt = tdr.estilo(interactivo=False)

    def test_ninguna_figura_toca_la_base(self):
        """Cada figura se dibuja con la carga de la base anulada: si alguna la
        necesitara, fallaría acá en vez de funcionar por accidente."""
        leer = pd.read_csv

        def solo_resultados(ruta, *a, **kw):
            if str(ruta).startswith(str(tdr.DB)):
                raise AssertionError(f"una figura leyó la base: {ruta}")
            return leer(ruta, *a, **kw)

        def prohibida(*a, **kw):
            raise AssertionError("una figura intentó cargar la base")

        dibujadas = 0
        with mock.patch.object(tdr, "cargar_db", prohibida), \
                mock.patch.object(tdr, "leer_aristas", prohibida), \
                mock.patch.object(pd, "read_csv", solo_resultados):
            for nombre, f in sorted(REG.items()):
                if f.faltan():
                    continue
                with self.subTest(figura=nombre):
                    fig = figuras.dibujar(nombre, guardar=False, plt=self.plt)
                    if fig is not None:
                        dibujadas += 1
                        self.plt.close(fig)
        if dibujadas == 0:
            self.skipTest("no hay tablas en resultados/: correr los notebooks")

    def test_una_figura_no_lee_lo_que_no_declaro(self):
        f = next(iter(REG.values()))
        t = figuras.Tablas(f)
        with self.assertRaises(KeyError):
            t("no_declarada.csv")


if __name__ == "__main__":
    unittest.main()
