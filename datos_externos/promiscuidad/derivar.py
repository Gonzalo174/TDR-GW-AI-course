#!/usr/bin/env python3
"""derivar.py — la lista de compuestos promiscuos, desde los datos crudos.

Es el único cálculo del proyecto que no se puede hacer desde `DB/`: contar
cuántas superestructuras contienen a cada compuesto exige las relaciones de
subestructura crudas, que la capa guardada en `DB/` ya no conserva (está
clusterizada). Por eso vive acá, fuera de los análisis, y su salida se versiona
como un insumo más, igual que `DB/`.

Necesita tres archivos privados (README.md de esta carpeta):

    TDR_RAW=<carpeta raw_data>  TDR_MAPEOS=<carpeta de mapeos>  python derivar.py

y escribe, en esta carpeta, las tres tablas que lee `analiceDB/03`:

    compuestos_promiscuos.csv   los que cumplen el criterio -> lo que filtra huerfanas/
    impacto_filtro.csv          cuántos compuestos y relaciones saca el filtro
    curva_filtrado.csv          sensibilidad del criterio (MW, umbral de promiscuidad)

No escribe la tabla por compuesto: el peso molecular con toda su precisión, al
lado del código, permitiría reidentificar los compuestos (PROVENANCE.md §3.13).
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

AQUI = Path(__file__).resolve().parent
RAIZ = AQUI.parents[1]
sys.path.insert(0, str(RAIZ / "comun"))
import tdr  # noqa: E402


def cargar_subestructuras_crudas(raw=None, mapeo=None, exigir_mapa=False):
    """Relaciones de subestructura crudas + pesos moleculares, mapeadas a `drug_id`.

    Lee `raw_data/subestructures_chembl35_biolip.txt` (1 065 346 relaciones) y
    `raw_data/compounds/compound_data.csv`. Es la unica entrada que no sale de
    `DB/`: la capa de subestructuras guardada ya esta clusterizada y no conserva
    el conteo de superestructuras por molecula.

    Los crudos traen los identificadores originales, mientras que `DB/` los tiene
    codificados. Sin traducir, cualquier cruce contra `datos.bioact` da vacio sin
    lanzar nada. Por eso se recodifican aca, con `mapa_compuesto`: se busca en
    `mapeo`, o en la variable de entorno TDR_MAPEOS. Si no esta, se avisa y se
    devuelven los ids originales, que solo sirven para mirar los crudos.

    Ni los crudos ni los mapeos se publican: en un clon del repositorio esta
    funcion no corre, y por eso su salida se versiona en esta carpeta.

    `exigir_mapa=True` convierte el aviso en error; es lo que usa `__main__`.
    Sin el mapa, la lista de promiscuos sale con 30 compuestos falsos en vez de
    los verdaderos, y nada falla (PROVENANCE.md §3.10).
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

    mc = _mapa_compuesto(mapeo)
    if mc is not None:
        for col in ("padre", "hijo"):
            sub[col] = sub[col].map(mc).fillna(tdr.FALTANTE).astype("int64")
        sub = sub[(sub["padre"] != tdr.FALTANTE) & (sub["hijo"] != tdr.FALTANTE)]
        comp["drug_id"] = comp["drug_id"].map(mc).fillna(tdr.FALTANTE).astype("int64")
        comp = comp[comp["drug_id"] != tdr.FALTANTE]
    elif exigir_mapa:
        raise FileNotFoundError(
            "falta mapa_compuesto.csv: definir TDR_MAPEOS con la carpeta de los mapeos "
            "privados. Sin el mapa los ids crudos no cruzan contra DB/ y la lista de "
            "promiscuos sale falsa (PROVENANCE.md, 3.10)")
    else:
        print("AVISO: sin mapa_compuesto; los drug_id quedan en la notacion "
              "original y NO cruzan contra DB/ (ver TDR_MAPEOS)", flush=True)
    return sub, comp


def _mapa_compuesto(mapeo=None):
    """{id original -> codigo} de `mapa_compuesto.csv`, o None si no esta."""
    import os
    d = Path(mapeo) if mapeo else Path(os.environ.get("TDR_MAPEOS", "/no-existe"))
    f = d / "mapa_compuesto.csv"
    if not f.exists():
        return None
    m = pd.read_csv(f)
    return dict(zip(m["original"].astype("int64"), m["codigo"].astype("int64")))


def promiscuidad_subestructural(sub, comp, mw_max=tdr.MW_PROMISCUIDAD,
                                n_par_min=tdr.N_PARENTALES_PROMISCUIDAD):
    """Compuestos promiscuos de la capa de subestructuras (S3 Fig del paper 2016).

    `N_parentales` es la cantidad de superestructuras que contienen a la
    molécula; el paper excluye las relaciones de las que además pesan menos de
    `mw_max` Da. Devuelve `(por_compuesto, promiscuos, impacto)`:

      por_compuesto  drug_id, n_parentales, molweight
      promiscuos     los que cumplen el criterio -> es lo que consume `huerfanas/`
      impacto        cuántos compuestos y cuántas aristas saca el filtro

    Es principal porque, aplicado (CONVENCIONES.md §10), cambia todos los números aguas
    abajo: semillas, cobertura y ranking de las pseudohuérfanas.
    """
    n_par = sub.groupby("hijo").size().rename("n_parentales")
    # compound_data trae filas repetidas para algunos compuestos: sin este
    # drop_duplicates el merge los duplica y la lista contaba 7 filas para 5
    # compuestos (en v4, 30 para 25). Corregido el 2026-09-24.
    mw = comp[["drug_id", "molweight"]].drop_duplicates("drug_id")
    d = (n_par.reset_index().rename(columns={"hijo": "drug_id"})
         .merge(mw, on="drug_id", how="left"))

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


if __name__ == "__main__":
    # el mapa se exige: sin el, los ids crudos no cruzan contra DB/ y la lista
    # sale con 30 compuestos falsos en vez de los verdaderos (PROVENANCE.md §3.10)
    sub, comp = cargar_subestructuras_crudas(exigir_mapa=True)
    por_compuesto, promiscuos, impacto = promiscuidad_subestructural(sub, comp)
    curva = curva_filtrado(por_compuesto, sub)
    promiscuos.to_csv(AQUI / "compuestos_promiscuos.csv", index=False)
    impacto.to_csv(AQUI / "impacto_filtro.csv", index=False)
    curva.to_csv(AQUI / "curva_filtrado.csv", index=False)
    print(f"{len(promiscuos)} compuestos promiscuos, sobre {len(sub)} relaciones -> {AQUI}")
