"""
funciones_analice.py — descriptivos de la base y de la red multicapa.

Casi todo auxiliar (conteos y figuras): el unico bloque principal son los
descriptivos que se citan como cifras en el paper —dimensiones de las tres capas,
disponibilidad por especie y el impacto del filtro de promiscuidad quimica—,
porque cambian numeros que van al texto.

Fuentes: `analiceDB/01_tablas_analisis.ipynb` y `02_conectividad_por_anotaciones.ipynb`
de v3, `reproducir_v6/descriptivos_red/` y `reproducir_v6/capa_quimica/02_filtro_promiscuidad.py`.

Nada se ejecuta al importar este modulo.
"""

from pathlib import Path

import numpy as np
import pandas as pd

import tdr


# ============================ PRINCIPALES ============================

def dimensiones_capas(datos, con_quimica=False):
    """Nodos y enlaces de las tres capas (Table 1-2 del paper 2019).

    Son las cifras que se citan en el texto, por eso es principal: un cambio en
    la carga de tablas se ve aca antes que en ningun otro lado.
    """
    f = [
        ("2 - proteínas", "proteínas (V_P)", datos.st["target_id"].nunique()),
        ("2 - proteínas", "especies", datos.st["sp_id"].nunique()),
        ("3 - afiliaciones", "categorías InterPro",
         datos.sta.loc[datos.sta["db"] == "ip", "ann"].nunique()),
        ("3 - afiliaciones", "categorías OrthoMCL",
         datos.sta.loc[datos.sta["db"] == "omcl", "ann"].nunique()),
        ("2-3 (afiliación)", "enlaces proteína-categoría", len(datos.sta)),
        ("1-2 (bioactividad)", "registros compuesto-proteína", len(datos.bioact)),
        ("1-2 (bioactividad)", "enlaces E_DP positivos", len(datos.posdt)),
        ("1-2 (bioactividad)", "enlaces E_DP negativos", len(datos.negdt)),
        ("1-2 (bioactividad)", "proteínas druggables", datos.posdt["target_id"].nunique()),
        ("1 - compuestos", "compuestos con bioactividad", datos.bioact["drug_id"].nunique()),
    ]
    if datos.bioact_org is not None:
        f.append(("fenotípica", "registros compuesto-organismo", len(datos.bioact_org)))
    if con_quimica and datos.tclus is not None:
        f += [
            ("1 - compuestos", "compuestos (V_D)", len(datos.tclus)),
            ("1 - compuestos", "clusters de fingerprint", datos.tclus["cluster_id"].nunique()),
            ("1 - compuestos", "clusters de subestructura", datos.sclus["cluster_id"].nunique()),
            ("1 (similitud)", "aristas E_DD tanimoto", len(datos.ddt)),
            ("1 (similitud)", "aristas E_DD subestructura", len(datos.dds)),
        ]
    return pd.DataFrame(f, columns=["capa", "entidad", "n"])


def disponibilidad_por_especie(datos):
    """Proteínas, anotaciones y druggables por especie (Table 3 del paper 2019)."""
    ip = datos.sta[datos.sta["db"] == "ip"]
    og = datos.sta[datos.sta["db"] == "omcl"]
    drg = set(datos.posdt["target_id"].astype(str))

    d = (datos.st.groupby("sp_id")["target_id"].nunique().rename("proteinas").reset_index()
         .merge(ip.groupby("sp_id")["target_id"].nunique().rename("con_interpro"),
                on="sp_id", how="left")
         .merge(og.groupby("sp_id")["target_id"].nunique().rename("con_orthomcl"),
                on="sp_id", how="left"))
    d["druggables"] = d["sp_id"].map(
        lambda s: len(set(datos.st.loc[datos.st["sp_id"] == s, "target_id"].astype(str)) & drg))
    d["pct_interpro"] = 100 * d["con_interpro"] / d["proteinas"]
    d["pct_orthomcl"] = 100 * d["con_orthomcl"] / d["proteinas"]
    d = d.merge(tdr.META, left_on="sp_id", right_on="sp", how="left").drop(columns=["sp"])
    return d.sort_values("druggables", ascending=False).reset_index(drop=True)


def grado_afiliaciones(datos):
    """Grado de cada categoría de afiliación (S1 Fig del paper 2016).

    El grado de una categoría es el número de proteínas anotadas a ella. La
    distribución es de cola pesada: unas pocas familias (quinasas,
    transportadores) concentran miles de proteínas, y son justamente las que la
    penalización por grado (`G'rk`, β = 1) atenúa.
    """
    g = (datos.sta.groupby(["ann", "db"])["target_id"].nunique()
         .reset_index(name="grado"))
    g["tramo"] = pd.cut(g["grado"], [0, 1, 5, 20, 100, 1000, np.inf],
                        labels=["1", "2-5", "6-20", "21-100", "101-1000", ">1000"])
    return g.sort_values("grado", ascending=False).reset_index(drop=True)


def promiscuidad_subestructural(sub, comp, mw_max=tdr.MW_PROMISCUIDAD,
                                n_par_min=tdr.N_PARENTALES_PROMISCUIDAD):
    """Compuestos promiscuos de la capa de subestructuras (S3 Fig del paper 2016).

    `N_parentales` es la cantidad de superestructuras que contienen a la
    molécula; el paper excluye las relaciones de las que además pesan menos de
    `mw_max` Da. Devuelve `(por_compuesto, promiscuos, impacto)`:

      por_compuesto  drug_id, n_parentales, molweight
      promiscuos     los que cumplen el criterio -> es lo que consume `huerfanas/`
      impacto        cuántos compuestos y cuántas aristas saca el filtro

    Es principal porque, aplicado (README §7.3), cambia todos los números aguas
    abajo: semillas, cobertura y ranking de las pseudohuérfanas.
    """
    n_par = sub.groupby("hijo").size().rename("n_parentales")
    d = (n_par.reset_index().rename(columns={"hijo": "drug_id"})
         .merge(comp[["drug_id", "molweight"]], on="drug_id", how="left"))

    sel = d[(d["molweight"] < mw_max) & (d["n_parentales"] > n_par_min)].copy()
    aristas = int(sub["hijo"].isin(set(sel["drug_id"])).sum())
    impacto = pd.DataFrame([{
        "criterio": f"MW < {mw_max} y N_parentales > {n_par_min}",
        "compuestos_filtrados": len(sel),
        "pct_compuestos": 100 * len(sel) / len(d),
        "aristas_filtradas": aristas,
        "pct_aristas": 100 * aristas / len(sub),
        "n_parentales_p99": float(d["n_parentales"].quantile(0.99)),
        "mw_mediana_filtrados": float(sel["molweight"].median()) if len(sel) else np.nan,
    }])
    return d, sel.sort_values("n_parentales", ascending=False).reset_index(drop=True), impacto


def composicion_bioactividades(datos):
    """Composición de las bioactividades por tag y por fuente.

    El tag `positive` es el que define la semilla de todo el modelo: su
    composición es la primera cifra que hay que poder defender.
    """
    cols = [c for c in ["activity_tag", "source", "cluster_consistent"]
            if c in datos.bioact.columns]
    return (datos.bioact.groupby(cols).size().reset_index(name="n")
            .assign(pct=lambda x: 100 * x["n"] / len(datos.bioact))
            .sort_values("n", ascending=False).reset_index(drop=True))


# ============================= AUXILIARES ============================

def tamanos_de_cluster(clusters):
    """Distribución de tamaños de cluster (fingerprint o subestructura)."""
    return (clusters.groupby("cluster_id").size().rename("tamano")
            .reset_index().sort_values("tamano", ascending=False).reset_index(drop=True))


def componentes_conexas(edges, con_bioactividad=None):
    """Componentes conexas de la capa química y cuáles tocan bioactividad positiva.

    `edges` son las aristas cluster-cluster (`clusID1`, `clusID2`);
    `con_bioactividad` el conjunto de clusters con alguna punta positiva.
    """
    import networkx as nx
    g = nx.Graph()
    g.add_edges_from(edges[["clusID1", "clusID2"]].values.tolist())
    comp = pd.DataFrame(
        {"CID": list(nx.connected_components(g)),
         "componente": np.arange(nx.number_connected_components(g))}
    ).explode("CID").reset_index(drop=True)
    comp["pos_bioact"] = comp["CID"].isin(con_bioactividad or set())
    return (comp.groupby("componente")
            .agg(total_clusters=("CID", "count"), clusters_con_bioact=("pos_bioact", "sum"))
            .reset_index())


def conectividad_entre_especies(datos, por="anotaciones"):
    """Matriz especie x especie: categorías (o drogas) compartidas.

    `por="anotaciones"` cuenta categorías de afiliación en común; `por="drogas"`
    cuenta compuestos con blanco positivo en las dos especies.
    """
    if por == "anotaciones":
        pares = datos.sta[["sp_id", "ann"]].drop_duplicates()
        clave = "ann"
    else:
        pares = (datos.posdt.merge(datos.st[["target_id", "sp_id"]], on="target_id")
                 [["sp_id", "drug_id"]].drop_duplicates())
        clave = "drug_id"
    sp = [s for s in tdr.ESPECIES_16 if s in set(pares["sp_id"])]
    conj = {s: set(pares.loc[pares["sp_id"] == s, clave]) for s in sp}
    m = pd.DataFrame([[len(conj[a] & conj[b]) for b in sp] for a in sp], index=sp, columns=sp)
    return m


def cargar_subestructuras_crudas(raw=None):
    """Relaciones de subestructura crudas + pesos moleculares, mapeadas a `drug_id`.

    Lee `raw_data/subestructures_chembl35_biolip.txt` (1 065 346 relaciones) y
    `raw_data/compounds/compound_data.csv`. Es la unica entrada de v4 que no sale
    de `DB/`: la capa de subestructuras guardada ya esta clusterizada y no
    conserva el conteo de superestructuras por molecula.
    """
    raw = Path(raw) if raw else tdr.RAW
    sub_p = raw / "subestructures_chembl35_biolip.txt"
    comp_p = raw / "compounds" / "compound_data.csv"
    if not sub_p.exists():
        raise FileNotFoundError(f"falta {sub_p}")

    sub = pd.read_csv(sub_p, sep="\t", low_memory=False)
    comp = pd.read_csv(comp_p, usecols=["compound_internal_id", "chembl_id", "molweight"],
                       dtype={"chembl_id": str}, low_memory=False)
    comp["drug_id"] = comp["compound_internal_id"].str.split("_").str[1].astype(int)

    if "substructure_internal_id" in sub.columns:            # ids internos
        sub["padre"] = sub["superstructure_internal_id"].str.split("_").str[1].astype("int64")
        sub["hijo"] = sub["substructure_internal_id"].str.split("_").str[1].astype("int64")
        if "same_compound" in sub.columns:
            sub = sub[~sub["same_compound"].astype(str).str.lower().eq("true")]
    else:                                                    # ids de ChEMBL
        cm = comp.dropna(subset=["chembl_id"]).drop_duplicates("chembl_id")
        mapa = dict(zip(cm["chembl_id"], cm["drug_id"]))
        sub["padre"] = sub["superstructure_id"].map(mapa)
        sub["hijo"] = sub["substructure_id"].map(mapa)

    sub = sub.dropna(subset=["padre", "hijo"])
    sub = sub[sub["padre"] != sub["hijo"]]
    sub[["padre", "hijo"]] = sub[["padre", "hijo"]].astype("int64")
    return sub, comp


def curva_filtrado(por_compuesto, sub, umbrales=(10, 50, 100, 500), mw_grilla=None):
    """Compuestos y aristas filtrables en función de (MW, umbral de promiscuidad).

    Sirve para ver si el criterio del paper cae en una zona estable o en el borde
    de un acantilado.
    """
    mw_grilla = mw_grilla if mw_grilla is not None else np.arange(50, 601, 25)
    filas = []
    for prom in umbrales:
        for mw in mw_grilla:
            sel = por_compuesto[(por_compuesto["molweight"] < mw) &
                                (por_compuesto["n_parentales"] > prom)]
            filas.append({"umbral_promiscuidad": prom, "mw_max": float(mw),
                          "compuestos_filtrables": len(sel),
                          "aristas_filtrables": int(sub["hijo"].isin(set(sel["drug_id"])).sum())})
    return pd.DataFrame(filas)


def escribir_meta(salidas, nb, **campos):
    """Reexporta `tdr.escribir_meta` (README §1.4)."""
    return tdr.escribir_meta(salidas, nb, **campos)


# --- figuras ---------------------------------------------------------------

def fig_disponibilidad(disp, plt):
    """Proteínas vs druggables por especie, con la recta de ajuste log-log."""
    d = disp[disp["sp_id"].isin(tdr.ESPECIES_16) & (disp["druggables"] > 0)]
    fig, ax = plt.subplots(figsize=(5.5, 4.5), tight_layout=True)
    ax.scatter(d["proteinas"], d["druggables"], color=tdr.S1, zorder=3)
    for _, r in d.iterrows():
        ax.annotate(tdr.NOMBRE_CORTO.get(r["sp_id"], r["sp_id"]),
                    (r["proteinas"] * 1.04, r["druggables"]), fontsize=7, color=tdr.INK2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("proteínas")
    ax.set_ylabel("druggables")
    ax.set_title("Disponibilidad de evidencia por especie")
    return fig


def fig_grado(grado, plt):
    """Distribución de grado de las categorías, por fuente (cola pesada)."""
    fig, ax = plt.subplots(figsize=(5.5, 4), tight_layout=True)
    bins = np.logspace(0, 4, 19)
#     bins = np.logspace(0, np.log10(max(grado["grado"].max(), 10)), 40)
    for db, etiqueta, color in [("ip", "InterPro", tdr.S1), ("omcl", "OrthoMCL", tdr.S2)]:
        g = grado.loc[grado["db"] == db, "grado"]
        ax.hist(g, bins=bins, histtype="step", color=color, label=f"{etiqueta} ({len(g)})")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("proteínas anotadas a la categoría")
    ax.set_ylabel("categorías")
    ax.legend(fontsize=8)
    return fig


def fig_conectividad(matriz, plt, titulo=""):
    """Matriz especie x especie de elementos compartidos."""
    fig, ax = plt.subplots(figsize=(6, 5), tight_layout=True)
    m = np.log10(matriz.values + 1)
    im = ax.imshow(m, cmap="Blues")
    etiquetas = [tdr.NOMBRE_CORTO.get(s, s) for s in matriz.columns]
    ax.set_xticks(range(len(etiquetas)))
    ax.set_xticklabels(etiquetas, rotation=90, fontsize=7)
    ax.set_yticks(range(len(etiquetas)))
    ax.set_yticklabels(etiquetas, fontsize=7)
    ax.grid(False)
    ax.set_title(titulo)
    
    # Crear el colorbar
    cbar = fig.colorbar(im, ax=ax, label="Cantidad compartida")
    
    # Configurar los ticks del colorbar manualmente
    manual_ticks = [0, 1, 10, 100, 1000, 10000]  # Valores originales deseados
    log_ticks = [np.log10(t + 1) for t in manual_ticks]  # Convertir a escala logarítmica
    cbar.set_ticks(log_ticks)  # Asignar los ticks al colorbar
    cbar.ax.set_yticklabels([f"{t}" for t in manual_ticks])  # Etiquetas manuales
    
    return fig


def fig_promiscuidad(por_compuesto, curva, plt):
    """MW vs N_parentales con el criterio del paper, y curva de filtrado."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), tight_layout=True)

    d = por_compuesto.dropna(subset=["molweight"])
    ax1.scatter(d["molweight"], d["n_parentales"], s=3, alpha=0.2, color=tdr.MUTED)
    sel = d[(d["molweight"] < tdr.MW_PROMISCUIDAD) &
            (d["n_parentales"] > tdr.N_PARENTALES_PROMISCUIDAD)]
    ax1.scatter(sel["molweight"], sel["n_parentales"], s=12, color=tdr.ST_CRIT,
                label=f"filtrados ({len(sel)})")
    ax1.axvline(tdr.MW_PROMISCUIDAD, color=tdr.S2, lw=1, ls="--")
    ax1.axhline(tdr.N_PARENTALES_PROMISCUIDAD, color=tdr.S2, lw=1, ls="--")
    ax1.set_yscale("log")
    ax1.set_xlabel("peso molecular (Da)")
    ax1.set_ylabel("superestructuras que lo contienen")
    ax1.set_xlim(0, 600)
    ax1.legend(fontsize=8)
    ax1.set_title("Criterio del paper (S3 Fig)")

    for umbral, color in zip(sorted(curva["umbral_promiscuidad"].unique()),
                             tdr.ORD + [tdr.S2]):
        c = curva[curva["umbral_promiscuidad"] == umbral]
        ax2.plot(c["mw_max"], 100 * c["aristas_filtrables"] / c["aristas_filtrables"].max(),
                 color=color, label=f"N_parentales > {umbral}")
    ax2.axvline(tdr.MW_PROMISCUIDAD, color=tdr.MUTED, lw=1, ls="--")
    ax2.set_xlabel("umbral de peso molecular (Da)")
    ax2.set_ylabel("aristas filtrables (% del máximo)")
    ax2.legend(fontsize=8)
    ax2.set_title("Sensibilidad al umbral")
    return fig
