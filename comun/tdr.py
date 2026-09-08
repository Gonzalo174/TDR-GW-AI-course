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

V4 = Path("/home/ggiordano/TDR/TDR_2026_v4")
DB = Path("/data1/TDR-2025/TDR_v7/gon3/DB")   # entrada; no se reescribe en v4
RAW = Path("/data1/TDR-2025/TDR_v7/raw_data")
GON4 = V4 / "gon4"                            # -> /data1/TDR-2025/TDR_v7/gon4

# Historico de v3: se conserva como control (README §1.3), ninguna notebook de
# v4 lo lee como insumo.
HISTORICO = Path("/data1/TDR-2025/TDR_v7/gon3/resultados")
CONTROL_V5 = HISTORICO / "genoma_completo_v5"

CARPETAS = ["createDB", "analiceDB", "genome_prioritization", "huerfanas"]

N_CORE_DEFAULT = 20

if str(V4 / "comun") not in sys.path:
    sys.path.insert(0, str(V4 / "comun"))

import nucleo as nf  # noqa: E402  (el modelo; ver comun/nucleo.py, congelado)


def out(carpeta):
    """Salidas del analisis `carpeta`: gon4/<carpeta>_out/ (la crea si no esta).

    Es el unico punto del arbol que escribe fuera de `gon4`, y no escribe nada:
    solo crea el directorio. Las salidas nunca van a `gon3/` ni al home.
    """
    if carpeta not in CARPETAS:
        raise ValueError(f"carpeta desconocida: {carpeta!r} (esperaba una de {CARPETAS})")
    d = GON4 / f"{carpeta}_out"
    (d / "figuras").mkdir(parents=True, exist_ok=True)
    return d


# ---------------------------------------------------------------------------
# Metadatos de especies
# ---------------------------------------------------------------------------

SPECIES = {
    "hsap": ("Homo sapiens",               "Mamíferos",     "Eucariota", False),
    "mmus": ("Mus musculus",               "Mamíferos",     "Eucariota", False),
    "atha": ("Arabidopsis thaliana",       "Plantas",       "Eucariota", False),
    "osat": ("Oryza sativa",               "Plantas",       "Eucariota", False),
    "scer": ("Saccharomyces cerevisiae",   "Hongos",        "Eucariota", False),
    "calb": ("Candida albicans",           "Hongos",        "Eucariota", True),
    "ecol": ("Escherichia coli",           "Bacterias",     "Procariota", True),
    "mtub": ("Mycobacterium tuberculosis", "Bacterias",     "Procariota", True),
    "sao":  ("Staphylococcus aureus",      "Bacterias",     "Procariota", True),
    "pfal": ("Plasmodium falciparum",      "Protozoos",     "Eucariota", True),
    "lmaj": ("Leishmania major",           "Protozoos",     "Eucariota", True),
    "tbrt": ("Trypanosoma brucei",         "Protozoos",     "Eucariota", True),
    "tcr":  ("Trypanosoma cruzi",          "Protozoos",     "Eucariota", True),
    "cele": ("Caenorhabditis elegans",     "Invertebrados", "Eucariota", False),
    "dmel": ("Drosophila melanogaster",    "Invertebrados", "Eucariota", False),
    "ddis": ("Dictyostelium discoideum",   "Amebozoos",     "Eucariota", False),
}
ESPECIES_16 = list(SPECIES)
KINETOPLASTIDOS = ["tcr", "tbrt", "lmaj"]

NOMBRE_CORTO = {
    "tcr": "T. cruzi", "tbrt": "T. brucei", "lmaj": "L. major", "pfal": "P. falciparum",
    "atha": "A. thaliana", "cele": "C. elegans", "ddis": "D. discoideum",
    "calb": "C. albicans", "osat": "O. sativa", "dmel": "D. melanogaster",
    "mmus": "M. musculus", "hsap": "H. sapiens", "sao": "S. aureus",
    "ecol": "E. coli", "mtub": "M. tuberculosis", "scer": "S. cerevisiae",
}

META = pd.DataFrame(
    [{"sp": k, "nombre": v[0], "grupo": v[1], "reino": v[2], "parasito": v[3]}
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
        a = a[(a["N"] > minimo) & a["sp_id"].notna() & (a["sp_id"] != "bioactive")]
        return a.sort_values("N", ascending=False)["sp_id"].tolist()


def cargar_db(anotaciones=True, quimica=False, fenotipo=False,
              cluster_consistent=True, verbose=True) -> Datos:
    """Carga las tablas de `DB/` (README: la base no se regenera en v4).

    anotaciones : capa 3 (InterPro tipo Domain + OrthoMCL) unida a la especie.
    quimica     : capa de compuestos (clusters + aristas). `01_edges_*` pesa 941 MB,
                  por eso esta apagada salvo pedido explicito.
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
    pos = d.bioact["activity_tag"] == "positive"
    neg = d.bioact["activity_tag"] == "negative"
    if cluster_consistent:
        pos &= d.bioact["cluster_consistent"]
        neg &= d.bioact["cluster_consistent"]
    d.posdt = d.bioact[pos].copy()
    d.negdt = d.bioact[neg].copy()

    if anotaciones:
        log("· 04a_interpro / 04b_orthomcl")
        ip = pd.read_csv(DB / "04a_interpro.csv")
        ip = ip[ip["type"] == "Domain"].drop_duplicates()
        og = pd.read_csv(DB / "04b_orthomcl.csv")
        d.sta = (pd.concat([ip[["target_id", "ann"]].assign(db="ip"),
                            og[og["ann"] != "-1"].assign(db="omcl")], ignore_index=True)
                 .merge(d.st[["target_id", "sp_id"]], on="target_id"))

    if quimica:
        log("· 01/02 clusters y aristas quimicas (pesado)")
        d.tclus = pd.read_csv(DB / "01_clusters_fingerprint.csv")
        d.sclus = pd.read_csv(DB / "02_clusters_subestructure.csv")
        d.ddt = pd.read_csv(DB / "01_edges_clusters_fingerprint.csv")
        d.dds = pd.read_csv(DB / "02_edges_clusters_subestructure.csv").set_index("from", drop=False)

    if fenotipo:
        log("· 03_bioactivities_organism_compound")
        d.bioact_org = pd.read_csv(DB / "03_bioactivities_organism_compound.csv")

    d.meta = fechas_db()
    return d


def fechas_db():
    """Fecha de modificacion de cada tabla de `DB/`, para el `meta.json` (§1.4)."""
    return {p.name: datetime.fromtimestamp(p.stat().st_mtime).strftime("%Y-%m-%d")
            for p in sorted(DB.glob("*.csv"))}


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

CACHE = V4 / "comun" / "datos_derivados"
PATRON_QUINASA = r"kinase|kinasa"


def nombres_ipr(refrescar=False):
    """Diccionario IPR -> descripcion, extraido de las salidas de InterProScan.

    Los .tsv de `raw_data/targets/` traen el accession IPR en la columna 12 y su
    descripcion en la 13. Se cachea en `comun/datos_derivados/ipr_nombres.csv`.
    """
    CACHE.mkdir(parents=True, exist_ok=True)
    cache = CACHE / "ipr_nombres.csv"
    if cache.exists() and not refrescar:
        d = pd.read_csv(cache)
        return dict(zip(d["ann"], d["nombre"]))

    import csv
    fuentes = list((RAW / "targets" / "interpro").glob("*.tsv"))
    fuentes += list((RAW / "targets" / "raw_target_data").glob("*.tsv"))
    nombres = {}
    for f in fuentes:
        with open(f, newline="", errors="replace") as fh:
            for row in csv.reader(fh, delimiter="\t"):
                if len(row) > 12 and row[11].startswith("IPR"):
                    nombres.setdefault(row[11], row[12])
    if not nombres:
        raise FileNotFoundError(f"sin .tsv de InterProScan bajo {RAW / 'targets'}")
    pd.DataFrame({"ann": list(nombres), "nombre": list(nombres.values())}).to_csv(cache, index=False)
    return nombres


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
