"""
tdr.py — lo compartido por las carpetas de analisis del repositorio.

Un unico lugar donde viven las rutas, la carga de las tablas, las metricas, la
paralelizacion y el estilo de las figuras (CONVENCIONES.md §1). La regla es: si una
funcion la necesitan dos carpetas, sube aca; nadie reimplementa una metrica.

  * rutas          -> `DB` (entrada, no se reescribe) y `out(carpeta)` (salidas
                      en `resultados/<carpeta>_out/`, CONVENCIONES.md §6)
  * carga          -> `cargar_db()` devuelve un `Datos` con las tablas pedidas
  * metricas       -> `pauc_normalizada`, `mcclish`, `auc01_mcclish`, `auc_global`
  * envoltorios    -> `semilla_global`, `relevance_scores`, `propagar`
  * paralelizacion -> `paralelizar(func, items, n_core=20)` con fork
  * figuras        -> `estilo()`, `guardar(fig, nombre, carpeta)`, paleta
  * procedencia    -> `escribir_meta(...)` deja el `NN_meta.json` de cada corrida

Uso desde la celda 1 de cualquier notebook, que busca la raiz del repositorio
(la carpeta que contiene `DB/`) y agrega `comun/` al path:

    RAIZ = Path.cwd()
    while not (RAIZ / "DB").is_dir() and RAIZ != RAIZ.parent:
        RAIZ = RAIZ.parent
    sys.path.insert(0, str(RAIZ / "comun"))
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
# Rutas (CONVENCIONES.md §6): cero literales de path dentro de los notebooks
# ---------------------------------------------------------------------------

# Todo cuelga de la raiz del repositorio, deducida de la ubicacion de este
# archivo: el repo se clona en cualquier lado y nada hay que reconfigurar.
RAIZ = Path(__file__).resolve().parents[1]
DB = RAIZ / "DB"                              # la base codificada; no se reescribe
SALIDAS = RAIZ / "resultados"                 # salidas de los notebooks

# Control independiente: los optimos por especie de la corrida v5, calculados
# con la metrica ya corregida (CONVENCIONES.md §3). Se versionan en el repo.
CONTROL_V5 = RAIZ / "control" / "genoma_completo_v5"

# Datos crudos (.tsv de InterProScan, tabla de subestructuras): NO se publican,
# son la base real, y ningun analisis los lee. Los usan `nombres_ipr` (que por
# eso no esta disponible) y `datos_externos/promiscuidad/derivar.py`, que
# calcula fuera de los analisis el unico insumo que no sale de DB/.
# Con una copia local se puede apuntar a ella: export TDR_RAW=/ruta/raw_data
RAW = Path(os.environ.get("TDR_RAW", RAIZ / "raw_data_ausente"))

CARPETAS = ["analiceDB", "genome_prioritization", "huerfanas"]

# Procesos para las corridas paralelas. Con TDR_NCORE se ajusta a la maquina:
# la corrida completa tarda ~45 min con 20 (ver README.md, "Tiempo de computo").
N_CORE_DEFAULT = int(os.environ.get("TDR_NCORE", 20))

if str(RAIZ / "comun") not in sys.path:
    sys.path.insert(0, str(RAIZ / "comun"))

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

# Etiquetas legibles de los tags, para ejes y tablas. Que un tag sea "positivo"
# no identifica nada: es la semantica del modelo, no notacion propia.
TAG_NOMBRE = {TAG_POSITIVE: "positivo", TAG_NEGATIVE: "negativo",
              TAG_INDETERMINATE: "indeterminado", TAG_INCONSISTENT: "inconsistente"}

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
# El codigo sigue siendo el identificador: es la clave con la que la especie
# entra en `DB/`, en `control/` y en las tablas de `resultados/`. Lo que se
# publica aca es a que organismo corresponde cada codigo, junto con el tipo de
# organismo, que el modelo usa —contrasta parasitos contra no parasitos, y
# procariotas contra eucariotas— y sin el cual las figuras por grupo no se leen.
#
# La correspondencia codigo -> organismo sale de los mapeos privados
# (`mapa_especie.csv` cruzado contra `raw_data/genomes/genome_data.csv`); esta
# tabla es la unica parte de esos diccionarios que se publica. El resto de la
# codificacion sigue en pie: compuestos, blancos, clusters, InterPro y OrthoMCL
# no se traducen, porque es lo que evita reidentificar los compuestos
# (PROVENANCE.md §3.13).

SPECIES = {
    # codigo: (organismo, grupo, reino, parasito)
     0: ("Mycobacterium tuberculosis", "Bacterias",     "Procariota", True),
     1: ("Escherichia coli",           "Bacterias",     "Procariota", True),
     3: ("Homo sapiens",               "Mamiferos",     "Eucariota",  False),
     4: ("Trypanosoma cruzi",          "Protozoos",     "Eucariota",  True),
     5: ("Staphylococcus aureus",      "Bacterias",     "Procariota", True),
     9: ("Caenorhabditis elegans",     "Invertebrados", "Eucariota",  False),
    10: ("Saccharomyces cerevisiae",   "Hongos",        "Eucariota",  False),
    15: ("Trypanosoma brucei",         "Protozoos",     "Eucariota",  True),
    17: ("Drosophila melanogaster",    "Invertebrados", "Eucariota",  False),
    18: ("Mus musculus",               "Mamiferos",     "Eucariota",  False),
    21: ("Arabidopsis thaliana",       "Plantas",       "Eucariota",  False),
    22: ("Dictyostelium discoideum",   "Amebozoos",     "Eucariota",  False),
    23: ("Leishmania major",           "Protozoos",     "Eucariota",  True),
    25: ("Candida albicans",           "Hongos",        "Eucariota",  True),
    26: ("Plasmodium falciparum",      "Protozoos",     "Eucariota",  True),
    28: ("Oryza sativa",               "Plantas",       "Eucariota",  False),
}
ESPECIES_16 = list(SPECIES)

# Nombre completo del organismo, por codigo.
NOMBRE_ESPECIE = {c: e for c, (e, _g, _r, _p) in SPECIES.items()}

# Kinetoplastidos: el clado de interes del proyecto, tres de los protozoos
# parasitos (T. cruzi, T. brucei, L. major). Se los trata como grupo en varias
# figuras.
KINETOPLASTIDOS = [4, 15, 23]

# La especie sobre la que se aplica el modelo en `huerfanas/04`: Plasmodium
# falciparum, agente de la malaria.
SP_FOCO = 26

# Etiqueta para ejes y tablas: el binomio abreviado ("P. falciparum"), que es
# lo que entra en un eje sin desbordarlo. El grupo y el reino viajan en columnas
# aparte de META, asi que no hacen falta en la etiqueta.
NOMBRE_CORTO = {c: f"{e.split()[0][0]}. {' '.join(e.split()[1:])}"
                for c, e in NOMBRE_ESPECIE.items()}

META = pd.DataFrame(
    [{"sp": k, "nombre": NOMBRE_CORTO[k], "grupo": v[1], "reino": v[2], "parasito": v[3]}
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
                  cluster. `huerfanas/` lo usa (como en la versión v3);
                  el barrido de `genome_prioritization/` no (como la corrida de control).
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
# Metricas — el unico lugar donde vive la pAUC (CONVENCIONES.md §1)
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
# Procedencia: el meta.json de cada corrida (CONVENCIONES.md §4)
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
    notebook: `f"{NB}_f01_<nombre>"` (CONVENCIONES.md §7)."""
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
    identifican por su codigo y no se pueden nombrar. La unica excepcion es
    `ANN_QUINASA`, donde se publica el accession de 391 dominios porque
    `genome_prioritization/03` los necesita (PROVENANCE.md §3.15).

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
    out = pd.concat(partes, ignore_index=True)
    # El sufijo sale del nombre del archivo y por lo tanto es texto. Cuando son
    # codigos de especie hay que devolverlos enteros: si no, `especie == SP_FOCO`
    # compara "26" contra 26, no encuentra nada y el notebook aborta.
    if out[columna].str.fullmatch(r"-?\d+").all():
        out[columna] = out[columna].astype("int64")
    return out


# ---------------------------------------------------------------------------
# Filtro de promiscuidad quimica (CONVENCIONES.md §10, decidido: se aplica)
# ---------------------------------------------------------------------------
#
# El paper 2016 (S3 Fig) excluye las relaciones de subestructura de las moleculas
# chicas y promiscuas: `MW < 150 Da` y `N_parentales > 100`, donde N_parentales es
# la cantidad de superestructuras que contienen a la molecula.
#
# La lista se calcula desde datos crudos en `datos_externos/promiscuidad/`;
# `analiceDB/03` la transcribe y mide su efecto sobre la base, y `huerfanas/` la
# aplica antes de construir cualquier semilla (CONVENCIONES.md §10).

MW_PROMISCUIDAD = 150       # Da
N_PARENTALES_PROMISCUIDAD = 100


def compuestos_promiscuos(salidas_analice=None):
    """Lee `analiceDB_out/03_compuestos_promiscuos.csv` y devuelve los `drug_id`.

    `analiceDB/03_calidad_bioactividades.ipynb` la transcribe desde
    `datos_externos/promiscuidad/`: hay que correrlo antes que `huerfanas/`. Si
    falta, el error dice exactamente eso.
    """
    salidas = Path(salidas_analice) if salidas_analice else (SALIDAS / "analiceDB_out")
    p = salidas / "03_compuestos_promiscuos.csv"
    if not p.exists():
        raise FileNotFoundError(
            f"falta {p}: correr analiceDB/03_calidad_bioactividades.ipynb antes "
            f"que huerfanas/ (CONVENCIONES.md §10, el filtro se aplica)")
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


# ---------------------------------------------------------------------------
# Quinasas: la categoria grande y promiscua contra la que se mide beta
# ---------------------------------------------------------------------------
#
# Cuales de los 14 083 dominios son de quinasa no se puede decidir desde `DB/`:
# hace falta el nombre, y `ann` es un entero codificado (ver `nombres_ipr`).
#
# Revelacion parcial y deliberada (PROVENANCE.md §3.15): se publica el
# accession de InterPro de estos 391 dominios al lado de su codigo. Es lo unico
# del diccionario privado `mapa_interpro.csv` que sale a la luz; el resto de los
# 36 924 accessions sigue sin publicarse. De aca para abajo el analisis usa
# **solo** el codigo: el accession esta para que se pueda auditar la marca, no
# para que el codigo lo use.
#
# Salen del catalogo de InterPro (`raw_data/targets/interpro.xml`): las entradas
# cuyo nombre o nombre corto contiene `kinase`, de tipo `Domain` —el mismo
# recorte que hace `cargar_db` con `IPR_DOMAIN`— y presentes en la base. Se
# regeneran con `datos_externos/quinasas/derivar.py`.

PATRON_QUINASA_FUENTE = r"kinase"   # el patron con el que se filtro el catalogo

ANN_QUINASA = {
    15: "IPR047895", 187: "IPR015795", 202: "IPR015865", 329: "IPR025287",
    407: "IPR041747", 415: "IPR035448", 609: "IPR047314", 651: "IPR034669",
    1024: "IPR054327", 1044: "IPR024637", 1290: "IPR033928", 1326: "IPR056301",
    1529: "IPR011147", 1546: "IPR057756", 1791: "IPR037310", 1829: "IPR018941",
    2254: "IPR015725", 2255: "IPR055214", 2336: "IPR032834", 2360: "IPR046854",
    2645: "IPR010622", 2793: "IPR059117", 2903: "IPR032640", 3024: "IPR049571",
    3084: "IPR035066", 3132: "IPR003117", 3192: "IPR027084", 3370: "IPR042709",
    3374: "IPR034848", 3387: "IPR015966", 3575: "IPR039154", 3598: "IPR020635",
    3900: "IPR034665", 3979: "IPR026611", 4051: "IPR022488", 4174: "IPR004399",
    4212: "IPR046855", 4273: "IPR000023", 4479: "IPR040712", 4566: "IPR030220",
    4573: "IPR035850", 4762: "IPR057764", 4914: "IPR019539", 4950: "IPR034667",
    4965: "IPR058710", 4984: "IPR015275", 5255: "IPR041087", 5300: "IPR035765",
    5312: "IPR008144", 5378: "IPR013727", 5530: "IPR045579", 5565: "IPR042704",
    5713: "IPR037317", 5812: "IPR007371", 5816: "IPR057529", 5818: "IPR037952",
    5844: "IPR054017", 6120: "IPR035771", 6349: "IPR002219", 6471: "IPR042666",
    6899: "IPR040999", 6960: "IPR031314", 7023: "IPR000961", 7144: "IPR056802",
    7319: "IPR029099", 7382: "IPR035063", 7456: "IPR033702", 7507: "IPR011611",
    7553: "IPR048002", 7761: "IPR004147", 7806: "IPR041108", 7809: "IPR042785",
    8001: "IPR004166", 8140: "IPR034668", 8183: "IPR047588", 8217: "IPR039137",
    8273: "IPR058126", 8335: "IPR045581", 8359: "IPR011712", 8428: "IPR029304",
    8596: "IPR012582", 8611: "IPR047368", 8719: "IPR035056", 8830: "IPR037778",
    8928: "IPR007862", 8954: "IPR030484", 9017: "IPR002420", 9036: "IPR035852",
    9122: "IPR055495", 9138: "IPR037709", 9164: "IPR041744", 9183: "IPR047908",
    9198: "IPR047256", 9245: "IPR044767", 9301: "IPR040665", 9303: "IPR042714",
    9505: "IPR034676", 9520: "IPR047916", 9692: "IPR047221", 9729: "IPR018984",
    9792: "IPR013750", 10000: "IPR039192", 10072: "IPR039482", 10131: "IPR034877",
    10169: "IPR034663", 10233: "IPR035870", 10349: "IPR035175", 10404: "IPR031634",
    10734: "IPR025198", 10778: "IPR041746", 11087: "IPR005702", 11152: "IPR024105",
    11196: "IPR025561", 11610: "IPR001048", 11667: "IPR054474", 11753: "IPR033768",
    11825: "IPR058209", 12060: "IPR042703", 12062: "IPR037708", 12081: "IPR035591",
    12232: "IPR034675", 12247: "IPR057469", 12252: "IPR013896", 12449: "IPR035805",
    12643: "IPR047484", 12684: "IPR037452", 12863: "IPR047475", 12868: "IPR047471",
    12920: "IPR037702", 12976: "IPR035589", 12992: "IPR003113", 13085: "IPR032380",
    13097: "IPR002192", 13259: "IPR037719", 13300: "IPR019380", 13339: "IPR022673",
    13597: "IPR037707", 13738: "IPR040976", 13777: "IPR037375", 13880: "IPR003594",
    13885: "IPR029597", 14037: "IPR001263", 14109: "IPR001245", 14300: "IPR000550",
    14392: "IPR013079", 14432: "IPR034673", 14453: "IPR018292", 14454: "IPR059233",
    14771: "IPR020676", 14879: "IPR000719", 14912: "IPR041739", 15002: "IPR017892",
    15092: "IPR057640", 15098: "IPR021820", 15261: "IPR001573", 15397: "IPR035022",
    15436: "IPR041734", 15445: "IPR044131", 15514: "IPR027916", 15546: "IPR004105",
    15622: "IPR057380", 15751: "IPR056782", 15819: "IPR049508", 15831: "IPR034662",
    15837: "IPR035064", 15965: "IPR035770", 16157: "IPR029878", 16236: "IPR034659",
    16358: "IPR024678", 16427: "IPR049871", 16469: "IPR034661", 16610: "IPR049761",
    16689: "IPR037705", 16696: "IPR014009", 16918: "IPR015285", 17091: "IPR047487",
    17100: "IPR034671", 17164: "IPR041745", 17195: "IPR030611", 17267: "IPR031994",
    17343: "IPR054000", 17376: "IPR001772", 17412: "IPR042743", 17477: "IPR041906",
    17617: "IPR044769", 17706: "IPR045270", 17901: "IPR022126", 17931: "IPR037706",
    17940: "IPR019247", 17985: "IPR015897", 18014: "IPR011102", 18155: "IPR024641",
    18265: "IPR003151", 18377: "IPR000403", 18438: "IPR045583", 18750: "IPR031775",
    18798: "IPR000341", 18893: "IPR031475", 19052: "IPR035574", 19088: "IPR035014",
    19100: "IPR047367", 19224: "IPR034677", 19293: "IPR049587", 19415: "IPR002826",
    19477: "IPR057614", 19545: "IPR010599", 19647: "IPR005467", 19661: "IPR014930",
    19686: "IPR040464", 19750: "IPR037606", 19953: "IPR021821", 20064: "IPR007373",
    20264: "IPR035692", 20463: "IPR042817", 20489: "IPR035583", 20603: "IPR006204",
    20626: "IPR037711", 20665: "IPR047222", 20676: "IPR037311", 20829: "IPR031831",
    20944: "IPR004358", 21092: "IPR023602", 21141: "IPR033470", 21184: "IPR042698",
    21363: "IPR029353", 21449: "IPR005189", 21475: "IPR011495", 21512: "IPR031782",
    21541: "IPR018485", 21542: "IPR042767", 21647: "IPR047965", 21683: "IPR031636",
    21775: "IPR058619", 21817: "IPR028754", 22053: "IPR047485", 22129: "IPR032872",
    22182: "IPR007521", 22204: "IPR047480", 22318: "IPR040867", 22409: "IPR012736",
    22483: "IPR041740", 22519: "IPR019510", 22623: "IPR022066", 22625: "IPR057292",
    22635: "IPR013543", 23019: "IPR026683", 23357: "IPR034907", 23650: "IPR033923",
    23666: "IPR035748", 23862: "IPR035804", 24126: "IPR035078", 24265: "IPR056383",
    24384: "IPR012844", 24404: "IPR040110", 24432: "IPR040642", 24439: "IPR035053",
    24465: "IPR057092", 24601: "IPR047499", 24674: "IPR042717", 24817: "IPR008145",
    24998: "IPR037784", 25034: "IPR041905", 25104: "IPR011126", 25130: "IPR035853",
    25136: "IPR042134", 25155: "IPR039430", 25401: "IPR045363", 25463: "IPR015793",
    25552: "IPR011104", 25646: "IPR042696", 25704: "IPR035751", 25787: "IPR040667",
    25951: "IPR047469", 26012: "IPR057564", 26050: "IPR002498", 26298: "IPR000687",
    26859: "IPR047896", 26980: "IPR042710", 27130: "IPR003175", 27149: "IPR057754",
    27272: "IPR032807", 27276: "IPR045495", 27497: "IPR034664", 27553: "IPR041390",
    27845: "IPR029462", 27864: "IPR035579", 27937: "IPR046861", 28115: "IPR016045",
    28153: "IPR048637", 28286: "IPR047486", 28299: "IPR022007", 28420: "IPR046803",
    28446: "IPR035020", 28489: "IPR024604", 28499: "IPR041309", 28586: "IPR031850",
    28587: "IPR057579", 28695: "IPR012737", 28843: "IPR034879", 29013: "IPR001206",
    29025: "IPR025200", 29036: "IPR041743", 29182: "IPR035572", 29342: "IPR035860",
    29377: "IPR035772", 29557: "IPR018955", 29647: "IPR039148", 29654: "IPR048629",
    29660: "IPR037704", 29713: "IPR022708", 29719: "IPR058681", 29796: "IPR003661",
    30038: "IPR035077", 30042: "IPR056392", 30087: "IPR024638", 30164: "IPR047477",
    30216: "IPR037716", 30270: "IPR042133", 30474: "IPR042718", 30493: "IPR054481",
    30850: "IPR022247", 30967: "IPR034670", 31036: "IPR019017", 31144: "IPR003852",
    31293: "IPR057465", 31519: "IPR006083", 31610: "IPR000756", 31655: "IPR045267",
    31660: "IPR056574", 31814: "IPR037638", 31862: "IPR013579", 31939: "IPR034674",
    32106: "IPR041328", 32198: "IPR035586", 32382: "IPR018484", 32390: "IPR015022",
    32535: "IPR035016", 32546: "IPR056865", 33007: "IPR044093", 33024: "IPR024585",
    33101: "IPR035851", 33113: "IPR029601", 33123: "IPR048470", 33227: "IPR042822",
    33279: "IPR035693", 33336: "IPR024953", 33340: "IPR011641", 33354: "IPR042697",
    33404: "IPR035060", 33422: "IPR047915", 33460: "IPR037313", 33495: "IPR022672",
    33555: "IPR034851", 33666: "IPR048394", 33741: "IPR042763", 33763: "IPR035837",
    33791: "IPR045067", 33792: "IPR042706", 33898: "IPR037312", 33907: "IPR034672",
    33924: "IPR019511", 33945: "IPR029477", 34087: "IPR028182", 34193: "IPR035588",
    34218: "IPR013695", 34239: "IPR010559", 34329: "IPR045801", 34557: "IPR054693",
    34640: "IPR008207", 34725: "IPR035533", 34756: "IPR033719", 34831: "IPR054466",
    34904: "IPR054521", 35003: "IPR047470", 35211: "IPR047962", 35249: "IPR039026",
    35259: "IPR047478", 35474: "IPR054352", 35487: "IPR056803", 35761: "IPR010737",
    35858: "IPR041429", 35873: "IPR013749", 35890: "IPR037703", 35896: "IPR044493",
    35983: "IPR042132", 36374: "IPR022049", 36401: "IPR056981", 36416: "IPR011620",
    36534: "IPR035062", 36787: "IPR011994", 36865: "IPR049870",
}


def anotaciones_quinasa():
    """Los `ann` de la capa de anotaciones que son dominios de quinasa.

    Vienen con el prefijo `IP` que `cargar_db` le repone a InterPro al
    concatenarlo con OrthoMCL (PROVENANCE.md §3.3), asi que cruzan directo
    contra `datos.sta` sin que el analisis tenga que saber de prefijos.
    """
    return np.array([f"IP{c}" for c in sorted(ANN_QUINASA)])


def marcar_quinasas(target_ids, sta, quinasas=None):
    """Serie booleana por proteina: True si tiene al menos un dominio de quinasa.

    `target_ids` fija el indice y el orden, asi que el resultado se puede pegar a
    un ranking sin reordenar nada. Las proteinas sin ninguna anotacion dan False.
    """
    q = set(map(str, quinasas if quinasas is not None else anotaciones_quinasa()))
    con_q = set(sta.loc[sta["ann"].astype(str).isin(q), "target_id"].astype(str))
    idx = pd.Index(map(str, target_ids), name="target_id")
    return pd.Series(idx.isin(con_q), index=idx, name="quinasa")
