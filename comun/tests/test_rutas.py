"""Convenciones de rutas y de numeracion (README §1.8 y §1.9).

Las dos reglas que sostienen el arbol de v4 y que es facil violar sin darse
cuenta: todas las salidas se resuelven contra el symlink `gon4` (nunca `gon3/`
ni el home), y todo archivo generado arranca con el numero del notebook que lo
escribio. Estos tests no crean directorios ni escriben nada en `gon4`.
"""
import sys, pathlib, unittest
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))   # comun/

import tdr as C


class TestRutas(unittest.TestCase):

    def test_gon4_apunta_a_data1(self):
        """El symlink `gon4` tiene que resolver a /data1: los resultados pesados
        no viven en el home."""
        self.assertTrue(C.GON4.is_symlink() or C.GON4.exists(),
                        f"falta el symlink {C.GON4}")
        self.assertTrue(str(C.GON4.resolve()).startswith("/data1"),
                        f"{C.GON4} resuelve a {C.GON4.resolve()}, no a /data1")

    def test_db_es_entrada_y_no_esta_bajo_gon4(self):
        """`DB/` es de solo lectura en v4: ninguna salida puede caer ahi."""
        self.assertNotIn(str(C.GON4.resolve()), str(C.DB))

    def test_out_solo_acepta_las_cuatro_carpetas(self):
        """Un typo en el nombre crearia una carpeta `_out` huerfana en gon4."""
        for mal in ("genome_priorization", "genome_prioritization/", "salidas", ""):
            with self.subTest(carpeta=mal):
                with self.assertRaises(ValueError):
                    C.out(mal)

    def test_las_carpetas_no_llevan_numero(self):
        """Un nombre que empieza con digito no es importable desde Python
        (README §1.1): el orden lo dan los numeros de los notebooks."""
        for c in C.CARPETAS:
            self.assertFalse(c[0].isdigit(), f"{c} arranca con digito")
            self.assertTrue(c.isidentifier(), f"{c} no es importable")

    def test_las_16_especies_estan_completas(self):
        self.assertEqual(len(C.ESPECIES_16), 16)
        self.assertEqual(set(C.ESPECIES_16), set(C.NOMBRE_CORTO))
        self.assertEqual(set(C.ESPECIES_16), set(C.META["sp"]))
        self.assertTrue(set(C.KINETOPLASTIDOS) <= set(C.ESPECIES_16))


class TestNumeracion(unittest.TestCase):
    """Los notebooks llevan numero y se lo propagan a todo lo que generan."""

    RAIZ = C.V4

    def test_los_notebooks_arrancan_con_numero(self):
        for carpeta in C.CARPETAS:
            for nb in sorted((self.RAIZ / carpeta).glob("*.ipynb")):
                with self.subTest(nb=nb.name):
                    self.assertTrue(nb.name[:2].isdigit(),
                                    f"{carpeta}/{nb.name} no arranca con NN_")

    def test_cada_carpeta_de_analisis_tiene_su_py(self):
        """Lo especifico de cada analisis vive en `funciones_<carpeta>.py`; lo
        compartido sube a `comun/` (README §1.1)."""
        for carpeta, py in [("analiceDB", "funciones_analice.py"),
                            ("genome_prioritization", "funciones_genome.py"),
                            ("huerfanas", "funciones_huerfanas.py")]:
            with self.subTest(carpeta=carpeta):
                self.assertTrue((self.RAIZ / carpeta / py).exists(),
                                f"falta {carpeta}/{py}")


if __name__ == "__main__":
    unittest.main()
