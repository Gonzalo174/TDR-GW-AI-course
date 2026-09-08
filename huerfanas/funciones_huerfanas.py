"""
funciones_huerfanas.py — pseudohuerfanas, cobertura de semilla y desorfanizacion.

Es `tdr-graph/orphan_drugs_v4.ipynb` mas `reproducir_v6/{cobertura_semilla,
desorfanizacion_global, inferencia_directa_indirecta, aplicacion_especie}` desarmados
en funciones. El ciclo que las combina vive en la celda de corrida de cada
notebook (README §1.2).

Division (README §5):
  PRINCIPAL = cambia un numero que va al paper. AUXILIAR = no.

**Filtro de promiscuidad**: por la decision del README §7.3 la capa quimica se
filtra ANTES de construir cualquier semilla, con `tdr.filtrar_capa_quimica()` y la
lista que deja `analiceDB/03`. Los numeros de esta carpeta no son comparables con
los de `gon3/resultados/pseudohuerfanas/`, que se calcularon sin filtrar.

Paralelizacion por fork: los workers leen las tablas de globals del modulo, que el
padre llena con `fijar_contexto()` antes de crear el pool.

Nada se ejecuta al importar este modulo.
"""

from pathlib import Path

import numpy as np
import pandas as pd

import tdr
import nucleo as nf

RSEED = 123457          # la semilla del muestreo de `orphan_drugs_v4.ipynb`
N_MUESTRA = 1000        # el paper 2016 toma 1 000 moleculas con un unico blanco

# Contexto compartido con los procesos hijos (se llena con `fijar_contexto`).
D = None            # Datos, con la capa quimica ya filtrada
POSDT_SP = None     # posdt con la especie del blanco anexada
CTX = {}


# ============================ PRINCIPALES ============================

def construir_pseudohuerfanas(datos, verbose=True):
    """Drogas pseudohuerfanas (k=1) y su estratificacion en grupos 0-3.

    Una droga es pseudohuerfana si tiene exactamente un blanco proteico conocido
    en una especie (`N == 1`, excluyendo `sp_id == "bioactive"`). Los grupos son

      0  sin enlaces quimicos a clusters con bioactividad positiva y blanco en
         una unica especie
      1  con al menos un enlace a un cluster con bioactividad positiva
      2  con blanco en mas de una especie
      3  ambas condiciones

    Requiere `cargar_db(quimica=True)`. Devuelve `(k1, posdt_sp)`, donde
    `posdt_sp` son los positivos con la especie del blanco anexada (se reutiliza
    para orfanizar en `priorizar_droga`).
    """
    if datos.tclus is None or datos.ddt is None:
        raise ValueError("hace falta cargar_db(quimica=True)")

    posdt_sp = datos.posdt.merge(datos.st[["target_id", "sp_id"]], on="target_id", how="left")
    # El merge `left` mete NaN y con eso `sp_id` pasa a float: los codigos de
    # especie tienen que volver a ser enteros o terminan en los nombres de
    # archivo como `01_26.0_pseudohuerfanas.csv`.
    posdt_sp = posdt_sp[posdt_sp["sp_id"].notna()].copy()
    posdt_sp["sp_id"] = posdt_sp["sp_id"].astype("int64")
    hits = (posdt_sp.groupby(["drug_id", "sp_id"]).target_id.count()
            .reset_index().rename(columns={"target_id": "N"})
            .merge(datos.tclus, left_on="drug_id", right_on="drug", how="left"))

    # clusters con alguna punta de bioactividad positiva
    e = datos.ddt[["clusID1", "clusID2"]].copy()
    con_bioact = set(hits["cluster_id"].dropna())
    e["p1"] = e["clusID1"].isin(con_bioact)
    e["p2"] = e["clusID2"].isin(con_bioact)
    alguno = e["p1"] | e["p2"]
    ee = pd.concat([e[alguno],
                    e[alguno].rename(columns={"clusID1": "clusID2", "clusID2": "clusID1",
                                              "p1": "p2", "p2": "p1"})])
    cid_edge_pos = set(ee.loc[ee["p2"], "clusID1"].unique())

    # Se excluye la pseudo-especie: agrupa blancos sin organismo asignado y no
    # es un genoma. Ojo al portar: aca decia `!= "bioactive"`, que contra los
    # codigos enteros de esta base es siempre verdadero y la dejaba entrar.
    k1 = hits[(hits["N"] == 1) & (hits["sp_id"] != tdr.SP_BIOACTIVE)].copy()
    n_esp = hits.groupby("cluster_id").sp_id.nunique()
    cid_multi = set(n_esp[n_esp > 1].index)

    k1["enlace_a_bioact"] = k1["cluster_id"].isin(cid_edge_pos)
    k1["multiespecie"] = k1["cluster_id"].isin(cid_multi)
    k1["grupo"] = 0
    k1.loc[k1["enlace_a_bioact"], "grupo"] = 1
    k1.loc[k1["multiespecie"], "grupo"] = 2
    k1.loc[k1["enlace_a_bioact"] & k1["multiespecie"], "grupo"] = 3

    if verbose:
        print(f"pseudohuerfanas: {len(k1)} drogas, {k1['sp_id'].nunique()} especies, "
              f"grupos {k1['grupo'].value_counts().sort_index().to_dict()}", flush=True)
    return k1.reset_index(drop=True), posdt_sp


def semilla_de_droga(datos, cid, posdt_sin_droga, tani_exp=1, sub_exp=1,
                     consolidar=False):
    """Semilla del vecindario quimico CSN(m) de una droga, pesada por similitud.

    Orfanizacion sin fuga: `posdt_sin_droga` ya no tiene NINGUNA arista de
    bioactividad de `cid`, asi que su blanco verdadero solo puede volver por otra
    droga vecina —que es la inferencia directa, y es legitima— nunca por si
    misma (ver `comun/tests/test_fuga.py`).

    `consolidar=False` devuelve la semilla sin consolidar por maximo, con la
    columna `type` (identity / tani / sub / direct): hace falta para clasificarla
    y para separar inferencia directa de indirecta. `consolidar=True` devuelve
    directamente lo que come `nucleo.nds()`.
    """
    drg = nf.get_druggable_targets(
        doi=[cid], posdt=posdt_sin_droga, st=datos.st, negdt=datos.negdt,
        use_neighbors=True, dds=datos.dds, ddt=datos.ddt,
        sclus=datos.sclus, tclus=datos.tclus,
        b_max_consolidate=consolidar, tani_exp=tani_exp, sub_exp=sub_exp)
    return drg.get(str(cid))


def clasificar_semilla(seed, blancos, sta):
    """nula / no informativa / informativa.

    - nula: no se encontro ningun vecino quimico con blanco conocido.
    - no informativa: hay semilla, pero ninguna de sus categorias funcionales
      (InterPro/OrthoMCL) coincide con las del blanco verdadero, con lo cual
      `nds()` le asigna score 0 de forma garantizada.
    - informativa: la semilla comparte al menos una categoria con el blanco.

    Es el cuello de botella del metodo: mientras el 61 % de las pseudohuerfanas
    tenga semilla nula, ninguna reponderacion mueve el numero.
    """
    if seed is None or len(seed) == 0:
        return "nula"
    cats_blanco = set(sta.loc[sta["target_id"].astype(str).isin(set(map(str, blancos))),
                              "ann"].dropna())
    cats_seed = set(sta.loc[sta["target_id"].isin(seed["target_id"].unique()),
                            "ann"].dropna())
    return "informativa" if (cats_blanco & cats_seed) else "no informativa"


def priorizar_droga(cid):
    """Una droga pseudohuerfana -> frank, rG, rSS, clase de semilla y de inferencia.

    Se ejecuta en un proceso hijo y lee `D`, `POSDT_SP` y `CTX` como globals.
    Una sola pasada por droga: del mismo `get_druggable_targets` sin consolidar
    salen la clasificacion de la semilla, la semilla consolidada para `nds()` y
    el desglose de evidencia que separa inferencia directa de indirecta.

      frank  posicion del blanco verdadero dentro de la especie, normalizada
      rG     posicion en el ranking global (todas las proteinas de la red)
      rSS    posicion dentro de la especie, sin normalizar
      clase  directa   -> el blanco verdadero esta entre las semillas
             indirecta -> solo se llega por la capa de afiliaciones
    """
    sp = CTX["sp_code"]
    filas_droga = POSDT_SP[POSDT_SP["drug_id"] == cid]
    blancos = filas_droga.loc[filas_droga["sp_id"] == sp, "target_id"].astype(str).unique()
    # orfanizacion: se remueven TODAS las aristas de bioactividad de la droga
    pposdt = POSDT_SP[POSDT_SP["drug_id"] != cid][["drug_id", "target_id"]]

    seed_raw = semilla_de_droga(D, cid, pposdt, tani_exp=CTX["tani_exp"],
                                sub_exp=CTX["sub_exp"], consolidar=False)
    semilla = clasificar_semilla(seed_raw, blancos, D.sta)

    base = {"drug_id": cid, "especie": sp, "semilla": semilla,
            "target_id": blancos[0] if len(blancos) else None,
            "n_semilla": 0 if seed_raw is None else seed_raw["target_id"].nunique()}

    if seed_raw is None or len(seed_raw) == 0:
        return {**base, "score": np.nan, "rank": np.nan, "frank": np.nan,
                "rG": np.nan, "rSS": np.nan, "clase": "sin semilla", "n_targets": np.nan}

    directa = len(set(seed_raw["target_id"].astype(str)) & set(blancos)) > 0
    seed = seed_raw.groupby("target_id", as_index=False)["w"].max()

    rnk = tdr.propagar(D.sta, seed, CTX["cat_rs"], beta=CTX["beta"],
                       lambda_=CTX["lambda_"], gamma=CTX["gamma"])
    rnk["target_id"] = rnk["target_id"].astype(str)
    rnk["rG"] = rnk["score"].rank(ascending=False, method="average")

    sp_rnk = rnk[rnk["target_id"].isin(CTX["sp_targets"])].copy()
    sp_rnk["rSS"] = sp_rnk["score"].rank(ascending=False, method="average")
    n_targets = len(sp_rnk)

    fila = sp_rnk[sp_rnk["target_id"].isin(blancos)]
    if len(fila) == 0:
        return {**base, "score": np.nan, "rank": np.nan, "frank": np.nan,
                "rG": np.nan, "rSS": np.nan,
                "clase": "directa" if directa else "indirecta", "n_targets": n_targets}

    f = fila.iloc[0]
    return {**base, "score": float(f["score"]), "rank": float(f["rSS"]),
            "frank": float(f["rSS"]) / n_targets,
            "rG": float(rnk.loc[rnk["target_id"] == f["target_id"], "rG"].iloc[0]),
            "rSS": float(f["rSS"]),
            "clase": "directa" if directa else "indirecta", "n_targets": n_targets}


def guardar_huerfanas(filas, salidas, nb, sp_code, sufijo="pseudohuerfanas"):
    """Escribe `<salidas>/<nb>_<especie>_<sufijo>.csv` (README §1.9)."""
    df = pd.DataFrame(filas)
    destino = Path(salidas) / f"{nb}_{sp_code}_{sufijo}.csv"
    df.to_csv(destino, index=False)
    if "frank" in df.columns and df["frank"].notna().any():
        print(f"  {sp_code}: {len(df)} drogas -> {destino.name}  "
              f"(frank<0.1 en {(df['frank'] < 0.1).mean():.1%}, "
              f"semilla {df['semilla'].value_counts().to_dict()})", flush=True)
    return df


def cobertura_por_palanca(datos, palanca, muestra, umbral=0.8, mapa_kegg=None,
                          verbose=True):
    """Cuánto sube la cobertura de semilla cada palanca del README §4.

    Las tres palancas sobre el 61 % de semilla nula, sin re-correr el modelo:
    lo que se mide es el techo que habilitaría cada una.

      "umbral"   umbral de similitud química. La capa guardada sólo conserva
                 pares con similitud >= 0.8: bajarlo requiere recalcular sobre
                 los fingerprints crudos, así que aquí sólo se puede SUBIR.
      "kegg"     KEGG como tercera fuente de afiliaciones (`mapa_kegg`: DataFrame
                 [target_id, ann]). Rescata semillas *no informativas*, no nulas:
                 agrega categorías, no vecinos.
      "fenotipo" la capa compuesto-organismo como semilla de respaldo. No da un
                 blanco proteico: da una especie, y la semilla se arma con sus
                 druggables. Es la información más débil que el paper descarta.

    Devuelve una fila por droga con la clase de semilla antes y después.
    """
    if palanca not in ("umbral", "kegg", "fenotipo"):
        raise ValueError(f"palanca desconocida: {palanca!r}")

    filas = []
    for _, r in muestra.iterrows():
        cid, sp = r["drug_id"], r["sp_id"]
        pposdt = datos.posdt[datos.posdt["drug_id"] != cid][["drug_id", "target_id"]]
        blancos = (datos.posdt[(datos.posdt["drug_id"] == cid)]["target_id"]
                   .astype(str).unique())

        d = datos
        if palanca == "umbral" and umbral > 0.8:
            import copy
            d = copy.copy(datos)
            d.ddt = datos.ddt[datos.ddt["weight"] >= umbral]

        seed = semilla_de_droga(d, cid, pposdt, consolidar=False)
        sta = datos.sta
        if palanca == "kegg" and mapa_kegg is not None:
            sta = pd.concat([datos.sta, mapa_kegg.assign(db="kegg")], ignore_index=True)
        antes = clasificar_semilla(seed, blancos, datos.sta)
        despues = clasificar_semilla(seed, blancos, sta)

        if palanca == "fenotipo" and antes == "nula" and datos.bioact_org is not None:
            # respaldo: los druggables de la especie contra la que el compuesto
            # tiene actividad fenotipica positiva
            org = datos.bioact_org
            # se seleccionan las bioactividades fenotípicas positivas
            sp_fen = org.loc[(org["drug_id"] == cid) &
                             (org["activity_tag"] == tdr.TAG_POSITIVE), "sp_id"].unique()
            if len(sp_fen):
                respaldo = pd.DataFrame({"target_id": np.concatenate(
                    [datos.druggables_de(s) for s in sp_fen]), "w": 1.0})
                despues = clasificar_semilla(respaldo, blancos, datos.sta)

        filas.append({"drug_id": cid, "especie": sp, "grupo": r.get("grupo"),
                      "palanca": palanca, "parametro": umbral if palanca == "umbral" else None,
                      "semilla_antes": antes, "semilla_despues": despues,
                      "rescatada": antes != "informativa" and despues == "informativa"})
    out = pd.DataFrame(filas)
    if verbose:
        print(f"palanca {palanca}: rescata {out['rescatada'].sum()} de {len(out)} "
              f"({out['rescatada'].mean():.1%})", flush=True)
    return out


def curva_recuperacion(rg, l_max=1000):
    """ρ(l) y λ(l) sobre las posiciones 1..l_max del ranking global (Fig 3A, 2016)."""
    rg = np.asarray(pd.to_numeric(pd.Series(rg), errors="coerce").dropna())
    l = np.arange(1, l_max + 1)
    rho = np.array([(rg <= x).sum() for x in l], dtype=float)
    lam = np.diff(rho, prepend=0.0)           # Δρ/Δl con Δl = 1
    return pd.DataFrame({"l": l, "rho": rho, "lambda": lam})


def r_g_estrella(curva, desde=50, k_sigma=3.0, ventana=5):
    """r*G: la posición donde λ(l) cae al ruido de fondo.

    λ(l) se suaviza con un spline cúbico (el paper usa una aproximación de tercer
    orden) y r*G es el primer l a partir de `desde` donde λ̃ baja del umbral
    λ∞ + k σ. Todo lo que quede por encima de r*G es la recuperación que el
    método atribuye a la señal y no al azar.
    """
    from scipy.interpolate import UnivariateSpline
    x = curva["l"].values.astype(float)
    y = (pd.Series(curva["lambda"].values.astype(float))
         .rolling(ventana, center=True, min_periods=1).mean().values)
    spl = UnivariateSpline(x, y, k=3, s=len(x) * np.var(y) * 0.05)
    suave = np.clip(spl(x), 0, None)

    cola = suave[x >= desde]
    lam_inf, sigma = float(cola.mean()), float(cola.std())
    umbral = lam_inf + k_sigma * sigma
    debajo = np.where((x >= desde) & (suave <= umbral))[0]
    rg_star = float(x[debajo[0]]) if len(debajo) else float(x[-1])

    curva = curva.copy()
    curva["lambda_suave"] = suave
    return curva, pd.DataFrame([{"lambda_inf": lam_inf, "sigma": sigma,
                                 "umbral": umbral, "r_g_estrella": rg_star,
                                 "k_sigma": k_sigma, "desde": desde}])


def embudo_especie(datos, sp=None):
    """Embudo de compuestos huérfanos con actividad fenotípica contra `sp`.

    `sp=None` usa `tdr.SP_FOCO`: el protozoo parásito sobre el que se aplica el
    modelo, identificado por su código.

      paso 1  compuestos con bioactividad fenotípica positiva contra `sp`
      paso 2  de esos, los que NO tienen ningún enlace de bioactividad a un
              blanco proteico (huérfanos en el sentido del paper)
      paso 3  de esos, los tratables: con al menos un vecino químico (Tanimoto o
              subestructura) que sí tenga blanco conocido

    Devuelve `(embudo, compuestos)`.
    """
    if datos.bioact_org is None:
        raise ValueError("hace falta cargar_db(fenotipo=True)")
    sp = tdr.SP_FOCO if sp is None else sp

    org = datos.bioact_org
    col_sp = "sp_id" if "sp_id" in org.columns else "organism_id"
    # se seleccionan las bioactividades fenotípicas positivas contra `sp`
    activos = set(org.loc[(org[col_sp] == sp) & (org["activity_tag"] == tdr.TAG_POSITIVE),
                          "drug_id"].unique())
    con_blanco = set(datos.bioact["drug_id"].unique())
    huerfanos = activos - con_blanco

    cid_con_blanco = set(datos.tclus.loc[datos.tclus["drug"].isin(con_blanco), "cluster_id"])
    tclus_h = datos.tclus[datos.tclus["drug"].isin(huerfanos)]
    vecinos = datos.ddt[datos.ddt["clusID1"].isin(cid_con_blanco) |
                        datos.ddt["clusID2"].isin(cid_con_blanco)]
    alcanzables = set(vecinos["clusID1"]) | set(vecinos["clusID2"])
    tratables = set(tclus_h.loc[tclus_h["cluster_id"].isin(alcanzables), "drug"])

    embudo = pd.DataFrame([
        {"paso": "1 - fenotípicamente activos", "n": len(activos)},
        {"paso": "2 - sin blanco proteico (huérfanos)", "n": len(huerfanos)},
        {"paso": "3 - tratables por el método", "n": len(tratables)},
    ])
    compuestos = pd.DataFrame({"drug_id": sorted(huerfanos)})
    compuestos["tratable"] = compuestos["drug_id"].isin(tratables)
    return embudo, compuestos


# ============================= AUXILIARES ============================

def muestra_pseudohuerfanas(k1, n=N_MUESTRA, rseed=RSEED, por_grupo=False):
    """Muestra aleatoria de drogas pseudohuerfanas.

    El paper 2016 toma 1 000 moleculas al azar con exactamente un blanco conocido
    (`por_grupo=False`). El pipeline v4 local muestrea hasta `n` por grupo, que es
    un protocolo distinto y mas exigente: se deja explicito para no confundir los
    numeros de una corrida con los del otro.
    """
    rng = np.random.default_rng(rseed)
    if not por_grupo:
        idx = rng.choice(len(k1), min(n, len(k1)), replace=False)
        return k1.iloc[np.sort(idx)].copy()
    partes = [g.iloc[np.sort(rng.choice(len(g), min(n, len(g)), replace=False))]
              for _, g in k1.groupby("grupo")]
    return pd.concat(partes, ignore_index=True)


def preparar_especie(datos, sp_code, params_opt, k1):
    """Contexto de una especie: RS con sus parametros optimos + sus targets.

    `params_opt` es `(alpha, beta, lambda_, gamma)` leido de
    `genome_prioritization_out/02_optimos_por_especie.csv`. A diferencia del
    barrido, aca la semilla NO es LOSO: cada droga se orfaniza por separado.
    """
    alpha, beta, lam, gamma = params_opt
    druggables = datos.posdt[["target_id"]].drop_duplicates().copy()
    druggables["w"] = 1
    cat_rs = tdr.relevance_scores(datos.sta, druggables, alpha=alpha)
    return {"sp_code": sp_code, "cat_rs": cat_rs[["ann", "RS"]],
            "beta": beta, "lambda_": lam, "gamma": gamma,
            "sp_targets": datos.st.loc[datos.st["sp_id"] == sp_code,
                                       "target_id"].astype(str).unique(),
            "tani_exp": 1, "sub_exp": 1,
            "n_pseudohuerfanas": int((k1["sp_id"] == sp_code).sum())}


def fijar_contexto(datos, posdt_sp, ctx):
    """Publica las tablas y el contexto como globals, para que los hijos del pool
    los vean por fork. Se llama en el padre, antes de `paralelizar`."""
    global D, POSDT_SP, CTX
    D, POSDT_SP, CTX = datos, posdt_sp, ctx


def paralelizar(func, items, n_core=tdr.N_CORE_DEFAULT, **kw):
    """Importada de `comun/tdr.py`."""
    return tdr.paralelizar(func, items, n_core=n_core, **kw)


def parametros_optimos(salidas_genome=None, nb="02"):
    """Optimos por especie calculados en `genome_prioritization/02`.

    Devuelve `{sp: (alpha, beta, lambda_, gamma)}`. Falla con un mensaje claro si
    el barrido todavia no se corrio: `huerfanas/` va despues (README §6.5).
    """
    salidas = Path(salidas_genome) if salidas_genome else (tdr.GON4 / "genome_prioritization_out")
    p = salidas / f"{nb}_optimos_por_especie.csv"
    if not p.exists():
        raise FileNotFoundError(
            f"falta {p}: correr genome_prioritization/01 y 02 antes que huerfanas/ "
            f"(README §6, pasos 3 y 5)")
    o = pd.read_csv(p)
    return {r["especie"]: (r["alpha"], r["beta"],
                           None if pd.isna(r["lambda_"]) else r["lambda_"], r["gamma"])
            for _, r in o.iterrows()}


def cargar_resultados(salidas, nb, patron="*"):
    """Concatena `<salidas>/<nb>_<patron>.csv`; la columna `especie` ya viene en
    las filas, asi que el sufijo del archivo se guarda como `archivo`."""
    return tdr.cargar_csvs(salidas, nb, patron=patron, columna="archivo")


def resumen_frank(df, corte=0.1):
    """% de pseudohuerfanas recuperadas (frank < corte), total, por grupo y por
    clase de semilla. Es la metrica con la que se compara cada iteracion."""
    def pct(d):
        return 100 * (d["frank"] < corte).mean()

    filas = [{"corte": corte, "particion": "total", "valor": "todas",
              "n": len(df), "pct_recuperadas": pct(df),
              "pct_sin_semilla": 100 * df["frank"].isna().mean()}]
    for col in ("grupo", "semilla", "clase", "especie"):
        if col not in df.columns:
            continue
        for v, d in df.groupby(col, dropna=False):
            filas.append({"corte": corte, "particion": col, "valor": v, "n": len(d),
                          "pct_recuperadas": pct(d),
                          "pct_sin_semilla": 100 * d["frank"].isna().mean()})
    return pd.DataFrame(filas)


def escribir_meta(salidas, nb, **campos):
    """Reexporta `tdr.escribir_meta` (README §1.4)."""
    return tdr.escribir_meta(salidas, nb, **campos)


# --- figuras ---------------------------------------------------------------

def fig_grupos(k1, plt):
    """Pseudohuérfanas por especie, apiladas por grupo 0-3."""
    t = (k1.groupby(["sp_id", "grupo"]).size().unstack(fill_value=0)
         .reindex(columns=[0, 1, 2, 3], fill_value=0))
    t = t.loc[t.sum(axis=1).sort_values(ascending=False).index]
    fig, ax = plt.subplots(figsize=(7, 3.5), tight_layout=True)
    abajo = np.zeros(len(t))
    etiquetas = {0: "ninguna", 1: "enlace a bioactividad +",
                 2: "multiespecie", 3: "ambas"}
    for g in [0, 1, 2, 3]:
        ax.bar(range(len(t)), t[g], bottom=abajo, color=tdr.COLOR_GRUPO[g],
               label=f"grupo {g} — {etiquetas[g]}")
        abajo += t[g].values
    ax.set_xticks(range(len(t)))
    ax.set_xticklabels([tdr.NOMBRE_CORTO.get(s, s) for s in t.index], rotation=90, fontsize=7)
    ax.set_ylabel("drogas pseudohuérfanas")
    ax.legend(fontsize=7)
    return fig


def fig_frank(df, plt, corte=0.1):
    """Distribución de frank por grupo (Fig 2 del paper 2016)."""
    fig, axes = plt.subplots(4, 1, figsize=(5.5, 6), sharex=True, sharey=True,
                             tight_layout=True)
    bins = np.arange(0, 0.61, 0.01)
    for g, ax in zip([0, 1, 2, 3], axes):
        d = df[df["grupo"] == g]
        ax.hist(d["frank"].dropna(), bins=bins, color=tdr.COLOR_GRUPO[g])
        ax.axvline(corte, color=tdr.INK2, lw=1, ls="--")
        ax.set_ylabel(f"grupo {g}", fontsize=8)
        ax.text(0.98, 0.85, f"{d['frank'].notna().sum()} de {len(d)} con semilla · "
                            f"{(d['frank'] < corte).mean():.0%} < {corte}",
                transform=ax.transAxes, ha="right", fontsize=7, color=tdr.INK2)
    axes[-1].set_xlabel("frank (posición relativa del blanco verdadero)")
    return fig


def fig_cobertura(cob, plt):
    """Clase de semilla por especie: el techo teórico del método."""
    t = (cob.groupby(["especie", "semilla"]).size().unstack(fill_value=0))
    for c in ("nula", "no informativa", "informativa"):
        if c not in t.columns:
            t[c] = 0
    t = t[["nula", "no informativa", "informativa"]]
    t = t.div(t.sum(axis=1), axis=0).sort_values("informativa")
    fig, ax = plt.subplots(figsize=(7, 3.5), tight_layout=True)
    abajo = np.zeros(len(t))
    for c in t.columns:
        ax.bar(range(len(t)), 100 * t[c], bottom=abajo, color=tdr.COLOR_SEMILLA[c], label=c)
        abajo += 100 * t[c].values
    ax.set_xticks(range(len(t)))
    ax.set_xticklabels([tdr.NOMBRE_CORTO.get(s, s) for s in t.index], rotation=90, fontsize=7)
    ax.set_ylabel("% de pseudohuérfanas")
    ax.set_title("Cobertura de semilla: informativa = techo del método")
    ax.legend(fontsize=8)
    return fig


def fig_palancas(pal, plt):
    """Ganancia potencial de cada palanca sobre la semilla nula / no informativa."""
    t = pal.groupby("palanca")["rescatada"].agg(["sum", "count", "mean"])
    fig, ax = plt.subplots(figsize=(5.5, 3.5), tight_layout=True)
    ax.bar(range(len(t)), 100 * t["mean"], color=tdr.CATEGORICOS[:len(t)])
    ax.set_xticks(range(len(t)))
    ax.set_xticklabels(t.index, fontsize=8)
    ax.set_ylabel("% de pseudohuérfanas rescatadas")
    for i, (n, c) in enumerate(zip(t["sum"], t["count"])):
        ax.text(i, 100 * n / c, f"{int(n)}/{int(c)}", ha="center", va="bottom", fontsize=7)
    return fig


def fig_recuperacion(curva, rg_star, plt):
    """Fig 3A: ρ(rG) y λ(rG) apilados, con la línea del umbral 3σ."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 5), sharex=True, tight_layout=True)
    ax1.plot(curva["l"], curva["rho"], color=tdr.S1)
    ax1.set_ylabel("ρ(l): drogas recuperadas")
    ax2.plot(curva["l"], curva["lambda"], color=tdr.MUTED, lw=1, label="λ(l)")
    ax2.plot(curva["l"], curva["lambda_suave"], color=tdr.S1, label="λ̃(l)")
    ax2.axhline(rg_star["umbral"].iloc[0], color=tdr.S2, ls="--", lw=1, label="λ∞ + 3σ")
    ax2.axvline(rg_star["r_g_estrella"].iloc[0], color=tdr.ST_CRIT, ls=":", lw=1.5,
                label=f"r*G = {rg_star['r_g_estrella'].iloc[0]:.0f}")
    ax2.set_xlabel("posición en el ranking global (rG)")
    ax2.set_ylabel("λ(l)")
    ax2.set_xscale("log")
    ax2.legend(fontsize=8)
    return fig


def fig_inferencia(df, plt):
    """Fig 4b: distribución de rG y rSS por clase de inferencia."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.5), tight_layout=True)
    clases = [c for c in ("directa", "indirecta") if c in set(df["clase"])]
    for ax, col, etiqueta in [(ax1, "rG", "rG (ranking global)"),
                              (ax2, "rSS", "rSS (dentro de la especie)")]:
        datos = [df.loc[df["clase"] == c, col].dropna() for c in clases]
        ax.boxplot(datos, showfliers=False)
        ax.set_xticks(range(1, len(clases) + 1))
        ax.set_xticklabels(clases, fontsize=8)
        ax.set_yscale("log")
        ax.set_ylabel(etiqueta)
    n = df["clase"].value_counts(normalize=True)
    ax1.set_title(" · ".join(f"{c}: {100*n.get(c, 0):.0f} %" for c in clases), fontsize=9)
    return fig


def fig_embudo(embudo, plt):
    """Compuestos que sobreviven a cada paso del embudo de la especie foco."""
    fig, ax = plt.subplots(figsize=(6, 3), tight_layout=True)
    ax.barh(range(len(embudo)), embudo["n"], color=tdr.ORD[:len(embudo)])
    ax.set_yticks(range(len(embudo)))
    ax.set_yticklabels(embudo["paso"], fontsize=8)
    ax.invert_yaxis()
    ax.set_xscale("log")
    for i, n in enumerate(embudo["n"]):
        ax.text(n * 1.05, i, f"{n:,}", va="center", fontsize=8, color=tdr.INK2)
    ax.set_xlabel("compuestos")
    return fig
