"""Integridad de las tablas de entrada (checklist del roadmap).

Estos tests leen `DB/` y tardan algunos minutos. Se saltean solos si la
base no está montada. Con `TDR_TESTS_QUIMICA=1` se agregan las verificaciones de
la capa química, que implican leer 941 MB de aristas.
"""
import os, sys, pathlib, unittest
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))   # comun/

import numpy as np
import pandas as pd

import tdr as C

HAY_DB = (C.DB / "00_specie_target.csv").exists()
CON_QUIMICA = os.environ.get("TDR_TESTS_QUIMICA") == "1"


@unittest.skipUnless(HAY_DB, "DB/ no disponible")
class TestIntegridad(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.d = C.cargar_db(anotaciones=True, quimica=CON_QUIMICA,
                            cluster_consistent=False, verbose=False)

    def test_proteinas_de_las_anotaciones_estan_en_la_tabla_de_especies(self):
        falt = set(self.d.sta["target_id"]) - set(self.d.st["target_id"])
        self.assertEqual(len(falt), 0, f"{len(falt)} proteínas anotadas sin especie")

    def test_druggables_estan_en_la_tabla_de_especies(self):
        falt = set(self.d.posdt["target_id"]) - set(self.d.st["target_id"])
        self.assertLess(len(falt) / self.d.posdt["target_id"].nunique(), 0.5,
                        "más de la mitad de los blancos druggables no mapean a especie")

    def test_prefijos_de_anotacion(self):
        """El Fisher de nds_fun separa por prefijo IP / OG: si aparece otro, los
        p-valores de esa categoría quedan en NaN sin aviso."""
        malos = self.d.sta.loc[~self.d.sta["ann"].str.startswith(("IP", "OG")), "ann"]
        self.assertEqual(len(malos.unique()), 0, f"prefijos inesperados: {malos.unique()[:5]}")

    def test_los_dos_vocabularios_no_se_pisan(self):
        """InterPro y OrthoMCL se codificaron por separado, cada uno desde 0: sin
        el prefijo que repone `cargar_db`, el mismo entero sería un dominio y un
        grupo de ortología a la vez y el Fisher los fusionaría."""
        ip = set(self.d.sta.loc[self.d.sta["db"] == "ip", "ann"])
        og = set(self.d.sta.loc[self.d.sta["db"] == "omcl", "ann"])
        self.assertEqual(len(ip & og), 0, f"{len(ip & og)} anotaciones compartidas")

    def test_coherencia_db_prefijo(self):
        ip = self.d.sta[self.d.sta["db"] == "ip"]
        og = self.d.sta[self.d.sta["db"] == "omcl"]
        self.assertTrue(ip["ann"].str.startswith("IP").all())
        self.assertTrue(og["ann"].str.startswith("OG").all())

    def test_sin_duplicados_en_la_capa_de_afiliacion(self):
        dup = self.d.sta.duplicated(subset=["target_id", "ann"]).sum()
        self.assertEqual(dup, 0, f"{dup} enlaces proteína-categoría duplicados")

    def test_tags_de_bioactividad_conocidos(self):
        """Los cuatro códigos de `activity_tag`, y ninguno más."""
        esperados = {C.TAG_POSITIVE, C.TAG_NEGATIVE,
                     C.TAG_INDETERMINATE, C.TAG_INCONSISTENT}
        self.assertTrue(set(self.d.bioact["activity_tag"].unique()) <= esperados)

    def test_el_filtro_de_positivos_no_quedo_vacio(self):
        """El filtro por código es la falla silenciosa del port: contra esta base
        `activity_tag == "positive"` no lanza nada y devuelve cero filas, y el
        pipeline entero sigue adelante produciendo NaN."""
        self.assertGreater(len(self.d.posdt), 0, "posdt vacío: revisá TAG_POSITIVE")
        self.assertGreater(len(self.d.negdt), 0, "negdt vacío: revisá TAG_NEGATIVE")

    def test_la_base_es_entera(self):
        """Invariante que declara `acondicionarDB`: ni una celda no entera."""
        for t in ("00_specie_target.csv", "03_bioactivities_target_compound.csv",
                  "04a_interpro.csv", "04b_orthomcl.csv"):
            df = pd.read_csv(C.DB / t, nrows=5000)
            noent = [c for c in df.columns if df[c].dtype.kind not in "iu"]
            self.assertEqual(noent, [], f"{t}: columnas no enteras {noent}")

    def test_positivos_y_negativos_disjuntos(self):
        pos = set(map(tuple, self.d.posdt[["drug_id", "target_id"]].values))
        neg = set(map(tuple, self.d.negdt[["drug_id", "target_id"]].values))
        self.assertEqual(len(pos & neg), 0, "hay pares droga-blanco positivos y negativos")

    def test_cada_proteina_tiene_una_sola_especie(self):
        n = self.d.st.groupby("target_id")["sp_id"].nunique()
        self.assertEqual(int((n > 1).sum()), 0)

    @unittest.skipUnless(CON_QUIMICA, "TDR_TESTS_QUIMICA != 1")
    def test_pesos_de_similitud_en_rango(self):
        w = self.d.ddt["weight"]
        self.assertGreaterEqual(w.min(), 0.0)
        self.assertLessEqual(w.max(), 1.0)
        self.assertGreaterEqual(w.min(), 0.8,
                                "la capa química debería estar filtrada a >= 0.8")

    @unittest.skipUnless(CON_QUIMICA, "TDR_TESTS_QUIMICA != 1")
    def test_cada_droga_en_un_solo_cluster(self):
        n = self.d.tclus.groupby("drug")["cluster_id"].nunique()
        self.assertEqual(int((n > 1).sum()), 0)


if __name__ == "__main__":
    unittest.main()
