"""Convenciones de rutas y de numeracion (CONVENCIONES.md §6 y §1.9).

Las dos reglas que sostienen el arbol y que es facil violar sin darse cuenta:
todas las rutas se resuelven contra la raiz del repositorio (ningun path
absoluto de la maquina donde se escribio el codigo), y todo archivo generado
arranca con el numero del notebook que lo escribio.
"""
import sys, pathlib, unittest
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))   # comun/

import tdr as C


class TestRutas(unittest.TestCase):

    def test_todo_cuelga_del_repositorio(self):
        """Nada de rutas absolutas de la maquina de origen: un clon tiene que
        funcionar sin reconfigurar nada."""
        for nombre in ("DB", "SALIDAS", "CONTROL_V5"):
            ruta = getattr(C, nombre)
            self.assertTrue(str(ruta).startswith(str(C.RAIZ)),
                            f"{nombre} = {ruta} cae fuera de {C.RAIZ}")

    def test_la_base_esta_y_esta_completa(self):
        """Las tablas de `DB/`, incluidas las dos partes de la capa de aristas."""
        self.assertTrue(C.DB.is_dir(), f"falta {C.DB}")
        for t in ("00_specie_target.csv", "03_bioactivities_target_compound.csv",
                  "04a_interpro.csv", "04b_orthomcl.csv", "meta.json"):
            self.assertTrue((C.DB / t).exists(), f"falta DB/{t}")
        partes = sorted(C.DB.glob("01_edges_clusters_fingerprint.part*.csv.gz"))
        self.assertEqual(len(partes), 2, f"esperaba 2 partes de aristas, hay {len(partes)}")

    def test_db_es_entrada_y_no_recibe_salidas(self):
        """`DB/` es de solo lectura: ninguna salida puede caer ahi."""
        self.assertNotIn(str(C.DB), str(C.SALIDAS))
        for c in ("analiceDB", "genome_prioritization", "huerfanas"):
            self.assertNotIn(str(C.DB), str(C.out(c)))

    def test_out_solo_acepta_las_carpetas_de_analisis(self):
        """Un error en el nombre crearia una carpeta `_out` huerfana en resultados/."""
        for mal in ("genome_priorization", "genome_prioritization/", "salidas", ""):
            with self.subTest(carpeta=mal):
                with self.assertRaises(ValueError):
                    C.out(mal)

    def test_las_carpetas_no_llevan_numero(self):
        """Un nombre que empieza con digito no es importable desde Python
        (CONVENCIONES.md §1): el orden lo dan los numeros de los notebooks."""
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

    RAIZ = C.RAIZ

    def test_los_notebooks_arrancan_con_numero(self):
        for carpeta in C.CARPETAS:
            for nb in sorted((self.RAIZ / carpeta).glob("*.ipynb")):
                with self.subTest(nb=nb.name):
                    self.assertTrue(nb.name[:2].isdigit(),
                                    f"{carpeta}/{nb.name} no arranca con NN_")

    def test_cada_carpeta_de_analisis_tiene_su_py(self):
        """Lo especifico de cada analisis vive en `funciones_<carpeta>.py`; lo
        compartido sube a `comun/` (CONVENCIONES.md §1)."""
        for carpeta, py in [("analiceDB", "funciones_analice.py"),
                            ("genome_prioritization", "funciones_genome.py"),
                            ("huerfanas", "funciones_huerfanas.py")]:
            with self.subTest(carpeta=carpeta):
                self.assertTrue((self.RAIZ / carpeta / py).exists(),
                                f"falta {carpeta}/{py}")


if __name__ == "__main__":
    unittest.main()
