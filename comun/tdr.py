"""
tdr.py — lo compartido por las cuatro carpetas de analisis de TDR_2026_v4.

Un unico lugar donde viven las rutas, la carga de las tablas, las metricas, la
paralelizacion y el estilo de las figuras (README §1.1). La regla es: si una
funcion la necesitan dos carpetas, sube aca; nadie reimplementa una metrica.

  * rutas          -> `DB` (entrada, no se reescribe) y `out(carpeta)` (salidas
                      en `gon4/<carpeta>_out/`, README §1.8)
  * carga          -> `cargar_db()` devuelve un `Datos` con las tablas pedidas
  * metricas       -> `pauc_normalizada`, `mcclish`, `auc01_mcclish`, `auc_global`
  * envoltorios    -> `semilla_global`, `relevance_scores`, `propagar`
  * paralelizacion -> `paralelizar(func, items, n_core=20)` con fork
  * figuras        -> `estilo()`, `guardar(fig, nombre, carpeta)`, paleta
  * procedencia    -> `escribir_meta(...)` deja el `NN_meta.json` de cada corrida

Uso desde la celda 1 de cualquier notebook:

    import sys
    V4 = "/home/ggiordano/TDR/TDR_2026_v4"
    sys.path.insert(0, f"{V4}/comun")
    import tdr, nucleo as nf

Nada se ejecuta al importar este modulo: `out()` es lo unico que crea directorios,
y solo cuando se lo llama.
"""

from __future__ import annotations

import json
import os
import sys
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Rutas (README §1.8): cero literales de path dentro de los notebooks
# ---------------------------------------------------------------------------

# Todo cuelga de la raiz del repositorio, deducida de la ubicacion de este
# archivo: el repo se clona en cualquier lado y nada hay que reconfigurar.
RAIZ = Path(__file__).resolve().parents[1]
V4 = RAIZ                                     # nombre historico, mismo directorio
DB = RAIZ / "DB"                              # la base codificada; no se reescribe
SALIDAS = RAIZ / "resultados"                 # salidas de los notebooks
GON4 = SALIDAS                                # nombre historico de `SALIDAS`

# Control independiente: los optimos por especie de la corrida v5, calculados
# con la metrica ya corregida (README §1.3). Se versionan en el repo.
CONTROL_V5 = RAIZ / "control" / "genoma_completo_v5"

# Datos crudos (.tsv de InterProScan, tabla de subestructuras): NO se publican,
# son la base real. Las dos funciones que los leen -`nombres_ipr` y
# `analiceDB.cargar_subestructuras_crudas`- quedan indisponibles en un clon.
# Con una copia local se puede apuntar a ella: export TDR_RAW=/ruta/raw_data
RAW = Path(os.environ.get("TDR_RAW", RAIZ / "raw_data_ausente"))

CARPETAS = ["createDB", "analiceDB", "genome_prioritization", "huerfanas"]

N_CORE_DEFAULT = 20

if str(V4 / "comun") not in sys.path:
    sys.path.insert(0, str(V4 / "comun"))

import nucleo as nf  # noqa: E402  (el modelo; ver comun/nucleo.py, congelado)


def out(carpeta):
    """Salidas del analisis `carpeta`: resultados/<carpeta>_out/ (la crea si no esta).

    Es el unico punto del arbol que decide donde se escribe, y no escribe nada:
    solo crea el directorio. Todo queda dentro del repositorio.
    """
    if carpeta not in CARPETAS:
        raise ValueError(f"carpeta desconocida: {carpeta!r} (esperaba una de {CARPETAS})")
    d = SALIDAS / f"{carpeta}_out"
    (d / "figuras").mkdir(parents=True, exist_ok=True)
    return d


# ---------------------------------------------------------------------------
# Metadatos de especies
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Codigos de la base
# ---------------------------------------------------------------------------
# La base de este repositorio esta codificada: toda columna es un entero, y los
# diccionarios que traducen esos enteros a la notacion original (UniProt,
# InterPro, OrthoMCL, ids de compuesto, nombres de especie) son privados y no se
# publican. Las constantes de abajo son las unicas equivalencias necesarias para
# leer el modelo, y estan escritas a mano justamente para no tener que publicar
# ningun diccionario.
#
# Cuidado al portar codigo de versiones anteriores: alli los filtros se escribian
# contra strings ("positive", "Domain"). Contra esta base esa comparacion no
# falla, devuelve cero filas en silencio.

# --- activity_tag ---
TAG_POSITIVE = 2       # se seleccionan las bioactividades positivas
TAG_NEGATIVE = 3       # bioactividades negativas
TAG_INDETERMINATE = 1  # medicion sin umbral de corte que la defina
TAG_INCONSISTENT = 0   # el mismo par compuesto-blanco medido con signos opuestos

# --- type, en 04a_interpro ---
IPR_DOMAIN = 5         # se selecciona el tipo Domain: la unidad funcional que
                       # usa el modelo. Los otros tipos de InterPro (Family,
                       # Homologous_superfamily, Active_site, Binding_site,
                       # Conserved_site, PTM, Repeat) no entran.

# --- origen, en 03_bioactivities_target_compound ---
ORIGEN_BIOLIP = 1
ORIGEN_CHEMBL = 2
ORIGEN_AMBOS = 0

# --- valores especiales ---
FALTANTE = -1          # todo faltante de la base, en cualquier columna
SP_BIOACTIVE = 14      # pseudo-especie: agrupa los blancos sin especie asignada
                       # que provienen del conjunto bioactivo. No es un organismo
                       # y queda fuera de toda cuenta por especie.
ESCALA_PESO = 100      # los pesos de las aristas van escalados x100 (0.82 -> 82)


# ---------------------------------------------------------------------------
# Las 16 especies del modelo
# ---------------------------------------------------------------------------
# Se identifican por su codigo, no por su nombre: la correspondencia
# codigo -> especie vive en los diccionarios privados. Lo que si se publica es
# el tipo de organismo, porque el modelo lo usa —contrasta parasitos contra no
# parasitos, y procariotas contra eucariotas— y porque sin el las figuras por
# grupo no se pueden leer.
#
# Advertencia de anonimato: tres de las 16 son la unica especie de su
# combinacion (grupo, parasito), asi que para ellas el par identifica al
# organismo. Es una consecuencia aceptada de publicar el grupo.

SPECIES = {
    # codigo: (grupo, reino, parasito)
     0: ("Bacterias",     "Procariota", True),
     1: ("Bacterias",     "Procariota", True),
     3: ("Mamiferos",     "Eucariota",  False),
     4: ("Protozoos",     "Eucariota",  True),
     5: ("Bacterias",     "Procariota", True),
     9: ("Invertebrados", "Eucariota",  False),
    10: ("Hongos",        "Eucariota",  False),
    15: ("Protozoos",     "Eucariota",  True),
    17: ("Invertebrados", "Eucariota",  False),
    18: ("Mamiferos",     "Eucariota",  False),
    21: ("Plantas",       "Eucariota",  False),
    22: ("Amebozoos",     "Eucariota",  False),
    23: ("Protozoos",     "Eucariota",  True),
    25: ("Hongos",        "Eucariota",  True),
    26: ("Protozoos",     "Eucariota",  True),
    28: ("Plantas",       "Eucariota",  False),
}
ESPECIES_16 = list(SPECIES)

# Kinetoplastidos: el clado de interes del proyecto, tres de los protozoos
# parasitos. Se los trata como grupo en varias figuras.
KINETOPLASTIDOS = [4, 15, 23]

# La especie sobre la que se aplica el modelo en `huerfanas/04`: un protozoo
# parasito, agente de una enfermedad desatendida. Se la nombra por su codigo.
SP_FOCO = 26

# Etiqueta para ejes y tablas: el codigo, y entre parentesis el tipo de
# organismo, que es lo que hace legible la figura.
NOMBRE_CORTO = {c: f"sp{c:02d} ({g}{', parasito' if par else ''})"
                for c, (g, _r, par) in SPECIES.items()}

META = pd.DataFrame(
    [{"sp": k, "nombre": NOMBRE_CORTO[k], "grupo": v[0], "reino": v[1], "parasito": v[2]}
     for k, v in SPECIES.items()])


# ---------------------------------------------------------------------------
# Carga de la base
# ---------------------------------------------------------------------------

@dataclass
class Datos:
    """Contenedor de las tablas del modelo. Los campos no pedidos quedan en None."""
    st: pd.DataFrame = None
    species: pd.DataFrame = None
    posdt: pd.DataFrame = None
    negdt: pd.DataFrame = None
    bioact: pd.DataFrame = None
    bioact_org: pd.DataFrame = None
    sta: pd.DataFrame = None
    tclus: pd.DataFrame = None
    sclus: pd.DataFrame = None
    ddt: pd.DataFrame = None
    dds: pd.DataFrame = None
    meta: dict = field(default_factory=dict)

    def targets_de(self, sp: str) -> np.ndarray:
        return self.st.loc[self.st["sp_id"] == sp, "target_id"].astype(str).unique()

    def druggables_de(self, sp: str) -> np.ndarray:
        pos = set(self.posdt["target_id"].astype(str))
        return np.array([t for t in self.targets_de(sp) if t in pos])

    def especies_con_druggables(self, minimo=10) -> list:
        """Especies con mas de `minimo` blancos druggables: las 16 del modelo."""
        a = self.st[self.st["target_id"].isin(self.posdt["target_id"].unique())]
        a = a.groupby("sp_id")["target_id"].nunique().reset_index(name="N")
        # se excluye la pseudo-especie: agrupa blancos sin organismo asignado
        a = a[(a["N"] > minimo) & a["sp_id"].notna() & (a["sp_id"] != SP_BIOACTIVE)]
        return a.sort_values("N", ascending=False)["sp_id"].tolist()


def cargar_db(anotaciones=True, quimica=False, fenotipo=False,
              cluster_consistent=True, escalar_peso=True, verbose=True) -> Datos:
    """Carga las tablas de `DB/` (README: la base no se regenera en v4).

    anotaciones : capa 3 (InterPro tipo Domain + OrthoMCL) unida a la especie.
    quimica     : capa de compuestos (clusters + aristas). `01_edges_*` son
                  36 152 622 aristas repartidas en dos `.gz`, por eso esta
                  apagada salvo pedido explicito. Ver `leer_aristas`.
    escalar_peso : devuelve `weight` en [0,1] como en las versiones anteriores.
                  En la base va como entero x100.
    fenotipo    : bioactividades compuesto-organismo.
    cluster_consistent : filtra los positivos/negativos inconsistentes a nivel
                  cluster. `huerfanas/` lo usa (como `orphan_drugs_v4.ipynb`);
                  el barrido de `genome_prioritization/` no (como full_genome_v5).
    """
    def log(*a):
        if verbose:
            print(*a, flush=True)

    d = Datos()
    log("· 00_specie_target")
    d.st = pd.read_csv(DB / "00_specie_target.csv")
    d.species = d.st[["sp_id", "taxon_id"]].drop_duplicates().reset_index(drop=True)

    log("· 03_bioactivities_target_compound")
    d.bioact = pd.read_csv(DB / "03_bioactivities_target_compound.csv")
    # se seleccionan las bioactividades positivas (y las negativas, que son el
    # contraste); `activity_tag` es un entero, ver TAG_* arriba
    pos = d.bioact["activity_tag"] == TAG_POSITIVE
    neg = d.bioact["activity_tag"] == TAG_NEGATIVE
    if cluster_consistent:
        # `cluster_consistent` es 0/1 en la base: se castea a bool para que el
        # `&` sea una conjuncion logica y no un and de bits
        cc = d.bioact["cluster_consistent"].astype(bool)
        pos &= cc
        neg &= cc
    d.posdt = d.bioact[pos].copy()
    d.negdt = d.bioact[neg].copy()

    if anotaciones:
        log("· 04a_interpro / 04b_orthomcl")
        ip = pd.read_csv(DB / "04a_interpro.csv")
        # se seleccionan las anotaciones de tipo Domain
        ip = ip[ip["type"] == IPR_DOMAIN].drop_duplicates()
        og = pd.read_csv(DB / "04b_orthomcl.csv")
        # Los dos vocabularios se codificaron por separado, cada uno arrancando
        # en 0: el entero 100 es a la vez un dominio InterPro y un grupo de
        # ortologia distintos (los 14 083 codigos de dominio caen dentro del
        # rango de OrthoMCL). Concatenarlos crudos fusionaria anotaciones que no
        # tienen nada que ver. El prefijo las vuelve a separar y ademas repone lo
        # que `nucleo.get_annot_druggability_pv` necesita: alli el test de Fisher
        # se corre por vocabulario y elige las filas con `ann.startswith("IP")`
        # y `("OG")`.
        ip = ip[["target_id", "ann"]].assign(ann=lambda x: "IP" + x["ann"].astype(str), db="ip")
        og = og[og["ann"] != FALTANTE]   # se descartan los blancos sin ortologia
        og = og[["target_id", "ann"]].assign(ann=lambda x: "OG" + x["ann"].astype(str), db="omcl")
        d.sta = (pd.concat([ip, og], ignore_index=True)
                 .merge(d.st[["target_id", "sp_id"]], on="target_id"))

    if quimica:
        log("· 01/02 clusters y aristas quimicas (pesado)")
        d.tclus = pd.read_csv(DB / "01_clusters_fingerprint.csv")
        d.sclus = pd.read_csv(DB / "02_clusters_subestructure.csv")
        d.ddt = leer_aristas(escalar=escalar_peso)
        d.dds = pd.read_csv(DB / "02_edges_clusters_subestructure.csv").set_index("from", drop=False)

    if fenotipo:
        log("· 03_bioactivities_organism_compound")
        d.bioact_org = pd.read_csv(DB / "03_bioactivities_organism_compound.csv")

    d.meta = fechas_db()
    return d


def leer_aristas(escalar=True):
    """La capa de aristas cluster-cluster de fingerprint, reconstruida.

    Va partida en `01_edges_clusters_fingerprint.part00.csv.gz` y `.part01...`
    porque GitHub rechaza archivos de mas de 100 MB pero no limita su cantidad:
    partirla publica las 36 152 622 aristas completas en lugar de recortarlas
    por umbral de peso. Leerla es concatenar los trozos.

    `escalar=True` devuelve `weight` en [0,1]; en la base va como entero x100.
    """
    partes = sorted(DB.glob("01_edges_clusters_fingerprint.part*.csv.gz"))
    if not partes:
        raise FileNotFoundError(f"sin partes de 01_edges_clusters_fingerprint en {DB}")
    e = pd.concat([pd.read_csv(x) for x in partes], ignore_index=True)
    if escalar:
        e["weight"] = e["weight"] / ESCALA_PESO
    return e


def fechas_db():
    """Fecha de modificacion de cada tabla de `DB/`, para el `meta.json` (§1.4)."""
    return {p.name: datetime.fromtimestamp(p.stat().st_mtime).strftime("%Y-%m-%d")
            for p in sorted(list(DB.glob("*.csv")) + list(DB.glob("*.csv.gz")))}


# ---------------------------------------------------------------------------
# Metricas — el unico lugar donde vive la pAUC (README §1.1)
# ---------------------------------------------------------------------------

FPR_MAX = 0.10


def pauc_normalizada(y_true, score, fpr_max=FPR_MAX):
    """pAUC(FPR<=fpr_max) dividida por fpr_max. Azar = fpr_max/2, perfecto = 1.

    Interpola en el borde: agrega el punto exacto (fpr_max, tpr(fpr_max)); truncar
    en el ultimo punto con fpr <= fpr_max subestima la pAUC.
    """
    from sklearn.metrics import roc_curve, auc as sk_auc
    y_true = np.asarray(y_true)
    if y_true.sum() == 0 or y_true.sum() == len(y_true):
        return float("nan")
    fpr, tpr, _ = roc_curve(y_true, score)
    tpr_max = np.interp(fpr_max, fpr, tpr)
    m = fpr <= fpr_max
    f = np.append(fpr[m], fpr_max)
    t = np.append(tpr[m], tpr_max)
    if len(f) < 2:
        return float("nan")
    return sk_auc(f, t) / fpr_max


def mcclish(pauc_norm, fpr_max=FPR_MAX):
    """Correccion de McClish sobre la pAUC ya normalizada por fpr_max.

    En esa escala el minimo (clasificador al azar) es fpr_max/2 y el maximo es 1:
        AUC01 = 1/2 * (1 + (p - pmin) / (pmax - pmin))
    de modo que azar -> 0.5 y clasificador perfecto -> 1.0.
    """
    p = np.asarray(pauc_norm, dtype=float)
    pmin, pmax = fpr_max / 2.0, 1.0
    return 0.5 * (1.0 + (p - pmin) / (pmax - pmin))


def auc01_mcclish(y_true, score, fpr_max=FPR_MAX):
    """AUC parcial corregida por McClish, calculada desde las etiquetas."""
    return mcclish(pauc_normalizada(y_true, score, fpr_max), fpr_max)


def auc_global(y_true, score):
    from sklearn.metrics import roc_auc_score
    y_true = np.asarray(y_true)
    if y_true.sum() in (0, len(y_true)):
        return float("nan")
    return roc_auc_score(y_true, score)


def youden(y_true, score):
    """Corte de Youden: max(TPR + (1 - FPR) - 1) sobre la ROC global."""
    from sklearn.metrics import roc_curve
    fpr, tpr, _ = roc_curve(y_true, score)
    return float(np.max(tpr + (1 - fpr) - 1))


def bootstrap_auc01(y_true, score, n_boot=2000, fpr_max=FPR_MAX, rseed=0):
    """Distribucion bootstrap de AUC01 (el paper reporta la significancia de las
    diferencias entre estrategias con 2 000 remuestreos)."""
    rng = np.random.default_rng(rseed)
    y_true = np.asarray(y_true)
    score = np.asarray(score)
    n = len(y_true)
    out_ = np.empty(n_boot)
    for b in range(n_boot):
        idx = rng.integers(0, n, n)
        out_[b] = auc01_mcclish(y_true[idx], score[idx], fpr_max)
    return out_


# ---------------------------------------------------------------------------
# Envoltorios sobre nucleo.py
# ---------------------------------------------------------------------------

def semilla_global(posdt, st, sp_out=None, peso_uniforme=True):
    """Vector semilla de druggables (todas las drogas), opcionalmente removiendo
    la evidencia de una especie (validacion leave-one-species-out)."""
    # `nucleo` envuelve `sp_out` en lista solo si es str, que era el caso cuando
    # las especies se nombraban con un codigo de texto. Ahora son enteros y hay que
    # envolverlos aca: un int suelto revienta en el `isin` de get_druggable_targets.
    if sp_out is not None and not isinstance(sp_out, (list, tuple, set, np.ndarray)):
        sp_out = [sp_out]
    drg = nf.get_druggable_targets(posdt=posdt, st=st, sp_out=sp_out)
    seed = list(drg.values())[0].copy()
    if peso_uniforme:
        seed["w"] = 1
    return seed


def relevance_scores(sta=None, seed=None, alpha=0.6, umbral=0.2, adjust=False, pv=None):
    """Fisher + relevance score.

    Si ya se calcularon los p-valores se pasan en `pv` y `sta`/`seed` sobran: el
    test de Fisher es la parte cara del pipeline y no depende de alpha, que solo
    reescala los p-valores ya calculados.
    """
    cat_rs = pv.copy() if pv is not None else nf.get_annot_druggability_pv(sta, seed)
    cat_rs["RS"] = nf.rs(cat_rs["pv"].values, adjust=adjust, umbral=umbral,
                         alpha=alpha, method="fdr_bh")
    return cat_rs


def propagar(sta, seed, cat_rs, beta=1.0, lambda_=1.0, gamma=1.0):
    """Propagacion (voting scheme). beta=0 => G'_r ; beta>0 => G'_rk."""
    lam = None if (lambda_ is None or (isinstance(lambda_, float) and np.isnan(lambda_))) else lambda_
    return nf.nds(sta, seed, lambda_=lam, hybrid_gamma=gamma,
                  cat_rs=cat_rs[["ann", "RS"]], beta=beta)


def ranking_especie(rnk, sp_targets, positivos):
    """Restringe el ranking global a una especie y marca los positivos."""
    r = rnk.copy()
    r["target_id"] = r["target_id"].astype(str)
    r = r[r["target_id"].isin(set(map(str, sp_targets)))].copy()
    r["score"] = r["score"].fillna(0.0)
    r["druggable"] = r["target_id"].isin(set(map(str, positivos))).astype(int)
    r["rank"] = r["score"].rank(ascending=False, method="average")
    return r.sort_values("score", ascending=False).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Paralelizacion
# ---------------------------------------------------------------------------

def paralelizar(func, items, n_core=N_CORE_DEFAULT, chunk=1, desc=None, verbose=True):
    """Aplica `func` a cada elemento de `items` en `n_core` procesos.

    Usa fork: las tablas grandes cargadas en el proceso padre quedan accesibles
    en los hijos sin volver a serializarlas. `func` debe ser una funcion de
    modulo (no un lambda) y leer esas tablas de variables globales.
    n_core=1 ejecuta en serie, util para depurar.
    """
    import concurrent.futures as cf
    import multiprocessing as mp

    items = list(items)
    if desc and verbose:
        print(f"[{desc}] {len(items)} tareas en {n_core} procesos", flush=True)
    if n_core <= 1:
        return [func(x) for x in items]

    ctx = mp.get_context("fork")
    res = []
    with cf.ProcessPoolExecutor(max_workers=n_core, mp_context=ctx) as ex:
        for i, r in enumerate(ex.map(func, items, chunksize=chunk)):
            res.append(r)
            if verbose and desc and (i + 1) % max(1, len(items) // 10) == 0:
                print(f"  {desc}: {i + 1}/{len(items)}", flush=True)
    return res


def un_hilo_por_proceso():
    """Fija a 1 los hilos de BLAS. Tiene que correr ANTES de importar numpy, por
    eso los notebooks lo repiten en la celda 1: aca queda solo como referencia."""
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
              "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ.setdefault(v, "1")


# ---------------------------------------------------------------------------
# Procedencia: el meta.json de cada corrida (README §1.4)
# ---------------------------------------------------------------------------

def escribir_meta(salidas, nb, notebook=None, params=None, n_core=None, **extra):
    """Escribe `<salidas>/<nb>_meta.json` con la procedencia de la corrida.

    Diez lineas de codigo que eliminan la ambiguedad que hoy obliga a
    reconstruir la historia a mano: fecha, notebook de origen, parametros,
    fecha de las tablas de `DB/` y `n_core`.
    """
    salidas = Path(salidas)
    meta = {
        "notebook": notebook or f"{nb}_?.ipynb",
        "nb": str(nb),
        "fecha": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "params": _jsonable(params) if params is not None else None,
        "n_core": n_core,
        "tablas_db": fechas_db(),
        "python": sys.version.split()[0],
        "pandas": pd.__version__,
    }
    meta.update({k: _jsonable(v) for k, v in extra.items()})
    destino = salidas / f"{nb}_meta.json"
    destino.write_text(json.dumps(meta, indent=2, ensure_ascii=False))
    print(f"  meta -> {destino}", flush=True)
    return meta


def leer_meta(salidas, nb):
    """Contrapartida de `escribir_meta`: devuelve el dict de `<nb>_meta.json`."""
    return json.loads((Path(salidas) / f"{nb}_meta.json").read_text())


def _jsonable(x):
    """Convierte a algo serializable (los None de la grilla de lambda incluidos)."""
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple, np.ndarray)):
        return [_jsonable(v) for v in x]
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, Path):
        return str(x)
    return x


# ---------------------------------------------------------------------------
# Figuras — paleta y estilo (instancia de referencia del sistema de dataviz)
# ---------------------------------------------------------------------------

SURFACE = "#fcfcfb"   # superficie del grafico
INK = "#0b0b0b"       # tinta primaria
INK2 = "#52514e"      # tinta secundaria
MUTED = "#898781"     # ejes y etiquetas
GRID = "#e1e0d9"      # grilla (hairline)
AXIS = "#c3c2b7"      # linea de base

# slots categoricos, en orden fijo (nunca ciclar)
S1, S2, S3, S4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
CATEGORICOS = [S1, S2, S3, S4]

# rampa ordinal (una sola tinta, claro -> oscuro)
ORD = ["#86b6ef", "#2a78d6", "#104281"]

# paleta de estado (reservada; siempre con icono + etiqueta, nunca color solo)
ST_GOOD, ST_WARN, ST_SERIOUS, ST_CRIT = "#0ca30c", "#fab219", "#ec835a", "#d03b3b"

# grupos de pseudohuerfanas (0-3), en el orden en que los define huerfanas/01
COLOR_GRUPO = {0: "#d03b3b", 1: "#eb6834", 2: "#eda100", 3: "#898781"}
# clases de semilla, de peor a mejor
COLOR_SEMILLA = {"nula": "#d03b3b", "no informativa": "#eda100", "informativa": "#1baf7a"}


def estilo(interactivo=True):
    """Configura matplotlib y devuelve el modulo pyplot.

    `interactivo=False` fuerza el backend Agg (scripts sin display); en notebook
    se deja el backend que ya haya elegido Jupyter.
    """
    import matplotlib
    if not interactivo:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "text.color": INK,
        "xtick.color": MUTED, "ytick.color": MUTED,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
        "axes.spines.top": False, "axes.spines.right": False,
        "font.size": 9, "axes.titlesize": 10, "axes.titleweight": "bold",
        "legend.frameon": False, "lines.linewidth": 2, "figure.dpi": 130,
    })
    return plt


def guardar(fig, nombre, carpeta, png=True):
    """Guarda pdf (+png) en `carpeta`. `nombre` arranca con el numero del
    notebook: `f"{NB}_f01_<nombre>"` (README §1.9)."""
    carpeta = Path(carpeta)
    carpeta.mkdir(parents=True, exist_ok=True)
    fig.savefig(carpeta / f"{nombre}.pdf", bbox_inches="tight")
    if png:
        fig.savefig(carpeta / f"{nombre}.png", bbox_inches="tight", dpi=200)
    print(f"  figura -> {carpeta / nombre}.pdf", flush=True)


# ---------------------------------------------------------------------------
# Nombres de anotaciones (para leer los rankings priorizados)
# ---------------------------------------------------------------------------

PATRON_QUINASA = r"kinase|kinasa"


def nombres_ipr(refrescar=False):
    """Diccionario de anotacion -> descripcion legible. NO disponible en el repo.

    Las descripciones salen de los .tsv de InterProScan de `raw_data/targets/`,
    que son la base real y no se publican, y vienen indexadas por el accession
    de InterPro (un accession `IPR……`), mientras que la columna `ann` de esta base es un
    entero codificado. Ligar uno con otro exige `mapa_interpro`, que es privado.

    Consecuencia: en un clon del repositorio los dominios de un resultado se
    identifican por su codigo y no se pueden nombrar. Es la decision abierta 1
    de PLAN.md.

    Con una copia local de los crudos y del mapeo se puede reconstruir:
    `export TDR_RAW=/ruta/raw_data` y pasar `mapeo` con el diccionario privado.
    """
    raise NotImplementedError(
        "nombres_ipr() necesita raw_data/ y mapa_interpro, que no se publican; "
        "ver la decision abierta 1 de PLAN.md")


def anotar_targets(target_ids, sta, ipr=None, max_ann=3):
    """Devuelve, por proteina, sus anotaciones y una descripcion legible."""
    ipr = ipr if ipr is not None else nombres_ipr()
    s = sta[sta["target_id"].astype(str).isin(set(map(str, target_ids)))].copy()
    s["nombre_ann"] = s["ann"].map(ipr).fillna(s["ann"])
    g = (s.groupby("target_id")["nombre_ann"]
          .apply(lambda x: "; ".join(pd.unique(x)[:max_ann])).reset_index(name="anotaciones"))
    g["es_quinasa"] = g["anotaciones"].str.contains(PATRON_QUINASA, case=False, regex=True)
    return g


# ---------------------------------------------------------------------------
# Lectura de salidas ya generadas (la celda 5 de los notebooks)
# ---------------------------------------------------------------------------

def cargar_csvs(salidas, nb, patron="*", columna="especie"):
    """Concatena `<salidas>/<nb>_<patron>.csv` agregando el sufijo como columna.

    Es lo que usa la celda 5 para leer lo que dejo la celda 4 sin recalcular
    nada: `cargar_csvs(SALIDAS, "01")` -> los 16 CSV del barrido con la especie.
    """
    salidas = Path(salidas)
    archivos = sorted(salidas.glob(f"{nb}_{patron}.csv"))
    if not archivos:
        raise FileNotFoundError(
            f"no hay {nb}_{patron}.csv en {salidas}: correr la celda de corrida primero")
    partes = [pd.read_csv(f).assign(**{columna: f.stem[len(nb) + 1:]}) for f in archivos]
    return pd.concat(partes, ignore_index=True)


# ---------------------------------------------------------------------------
# Filtro de promiscuidad quimica (README §7.3, decidido: se aplica)
# ---------------------------------------------------------------------------
#
# El paper 2016 (S3 Fig) excluye las relaciones de subestructura de las moleculas
# chicas y promiscuas: `MW < 150 Da` y `N_parentales > 100`, donde N_parentales es
# la cantidad de superestructuras que contienen a la molecula. La base local no lo
# aplica: 30 compuestos concentran 104 075 relaciones, el 12 % de la capa.
#
# `analiceDB/03` mide el efecto y escribe la lista; `huerfanas/` la aplica antes
# de construir cualquier semilla. Vive aca porque lo usan las dos carpetas (§1.1).

MW_PROMISCUIDAD = 150       # Da
N_PARENTALES_PROMISCUIDAD = 100


def compuestos_promiscuos(salidas_analice=None):
    """Lee `analiceDB_out/03_compuestos_promiscuos.csv` y devuelve los `drug_id`.

    Es la salida del notebook `analiceDB/03_calidad_bioactividades.ipynb`: hay que
    correrlo antes que `huerfanas/`. Si falta, el error dice exactamente eso.
    """
    salidas = Path(salidas_analice) if salidas_analice else (GON4 / "analiceDB_out")
    p = salidas / "03_compuestos_promiscuos.csv"
    if not p.exists():
        raise FileNotFoundError(
            f"falta {p}: correr analiceDB/03_calidad_bioactividades.ipynb antes "
            f"que huerfanas/ (README §7.3, el filtro se aplica)")
    return pd.read_csv(p)["drug_id"].astype("int64").unique()


def filtrar_capa_quimica(datos, promiscuos, verbose=True):
    """Saca a los compuestos promiscuos de la capa de subestructuras.

    La capa esta indexada por cluster, no por compuesto, asi que el filtro se
    aplica en dos pasos: se quitan esas drogas de `sclus` —con lo que dejan de
    ser alcanzables como vecinos al armar la semilla— y se quitan de `dds` las
    aristas cuyos clusters se quedaron sin ninguna droga. La capa de fingerprint
    (`tclus`/`ddt`) no se toca: el criterio del paper es sobre subestructuras.

    Modifica una copia: devuelve un `Datos` nuevo, `datos` queda intacto.
    """
    import copy
    if datos.sclus is None or datos.dds is None:
        raise ValueError("hace falta cargar_db(quimica=True)")

    d = copy.copy(datos)
    prom = set(np.asarray(promiscuos).tolist())
    n_sclus0, n_dds0 = len(d.sclus), len(d.dds)

    d.sclus = d.sclus[~d.sclus["drug"].isin(prom)].copy()
    vivos = set(d.sclus["cluster_id"].unique())
    d.dds = d.dds[d.dds["from"].isin(vivos) & d.dds["to"].isin(vivos)].copy()

    if verbose:
        print(f"filtro de promiscuidad (MW < {MW_PROMISCUIDAD} y "
              f"N_parentales > {N_PARENTALES_PROMISCUIDAD}): "
              f"{len(prom)} compuestos", flush=True)
        print(f"  sclus {n_sclus0} -> {len(d.sclus)} "
              f"({100 * (n_sclus0 - len(d.sclus)) / n_sclus0:.2f} % removido)", flush=True)
        print(f"  dds   {n_dds0} -> {len(d.dds)} "
              f"({100 * (n_dds0 - len(d.dds)) / n_dds0:.2f} % removido)", flush=True)
    return d
