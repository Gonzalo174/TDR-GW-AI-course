#!/usr/bin/env python3
"""derivar.py — de donde sale `tdr.ANN_QUINASA`.

Decidir si un dominio es de quinasa exige su **nombre**, y la columna `ann` de
`DB/` es un entero codificado: ligar uno con otro necesita `mapa_interpro.csv`,
que es privado. Es la limitacion que `tdr.nombres_ipr()` documenta.

La lista resultante **no** es un insumo versionado aparte: vive en
`comun/tdr.py` como `ANN_QUINASA`, con el accession de InterPro al lado de cada
codigo (PROVENANCE.md §3.15). Este script es la procedencia de esa constante:
la regenera y la imprime lista para pegar, para que la marca se pueda auditar y
rehacer sin que ningun notebook dependa de los datos crudos.

    TDR_RAW=<carpeta raw_data>  TDR_MAPEOS=<carpeta de mapeos>  python derivar.py

Imprime el bloque de `ANN_QUINASA` y los numeros de cobertura medidos contra
`DB/`. No escribe ningun archivo.
"""
import os
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

AQUI = Path(__file__).resolve().parent
RAIZ = AQUI.parents[1]
sys.path.insert(0, str(RAIZ / "comun"))
import tdr  # noqa: E402


def entradas_interpro(raw=None):
    """`interpro.xml` -> DataFrame (acc, short_name, type, name).

    Es el catálogo completo de InterPro (49 674 entradas), la fuente autoritativa
    del nombre y del tipo. Se prefiere a `protein2ipr.dat`, que trae lo mismo
    repetido por proteína y pesa 83 GB.
    """
    d = Path(raw) if raw else Path(os.environ.get("TDR_RAW", "/no-existe"))
    f = d / "targets" / "interpro.xml"
    if not f.exists():
        raise FileNotFoundError(
            f"falta {f}: definir TDR_RAW con la carpeta de los datos crudos. "
            f"Sin el catalogo no hay nombres de dominio y no se puede decidir "
            f"cual es quinasa (PROVENANCE.md §3.6)")

    filas = []
    for _ev, el in ET.iterparse(f, events=("end",)):
        if el.tag == "interpro":
            filas.append((el.get("id"), el.get("short_name"), el.get("type"),
                          el.findtext("name") or ""))
            el.clear()
    return pd.DataFrame(filas, columns=["acc", "short_name", "type", "name"])


def _mapa_interpro(mapeo=None):
    """{accession IPR… -> codigo de DB/} de `mapa_interpro.csv`."""
    d = Path(mapeo) if mapeo else Path(os.environ.get("TDR_MAPEOS", "/no-existe"))
    f = d / "mapa_interpro.csv"
    if not f.exists():
        raise FileNotFoundError(
            f"falta {f}: definir TDR_MAPEOS con la carpeta de los mapeos privados. "
            f"Sin el mapa los accessions no cruzan contra la columna `ann` de DB/")
    m = pd.read_csv(f)
    return dict(zip(m["original"], m["codigo"].astype("int64")))


def anotaciones_quinasa(entradas, mapa, patron=tdr.PATRON_QUINASA_FUENTE):
    """{codigo de `DB/` -> accession de InterPro} de los dominios de quinasa.

    Tres filtros, en este orden:

      1. el nombre o el nombre corto de InterPro contiene `kinase`;
      2. el tipo es `Domain`. Es el mismo recorte que hace la base
         (`tdr.IPR_DOMAIN`): las familias y las superfamilias homologas no entran
         en la capa de anotaciones, asi que marcarlas no cambiaria nada;
      3. el accession esta en el mapeo, o sea que la anotacion existe en `DB/`.

    Devuelve el par completo —codigo y accession— porque es lo que se pega en
    `tdr.ANN_QUINASA`: el codigo es lo que usa el analisis, el accession es lo
    que permite auditar la marca (PROVENANCE.md §3.15).
    """
    es_quinasa = (entradas["name"].str.contains(patron, case=False, regex=True, na=False) |
                  entradas["short_name"].str.contains(patron, case=False, regex=True, na=False))
    dom = entradas[es_quinasa & (entradas["type"] == "Domain")].copy()
    dom["codigo"] = dom["acc"].map(mapa)
    dom = dom[dom["codigo"].notna()]
    return {int(c): a for c, a in zip(dom["codigo"], dom["acc"])}


def impacto(quinasas, raiz=None):
    """Cuanto de la capa de anotaciones y del genoma alcanza la marca.

    Se mide contra `DB/` para que el numero quede en el repositorio aunque los
    crudos no esten: cuantas anotaciones de quinasa hay, a cuantas proteinas
    tocan y que fraccion del total representan.
    """
    raiz = Path(raiz) if raiz else tdr.RAIZ
    ip = pd.read_csv(raiz / "DB" / "04a_interpro.csv")
    ip = ip[ip["type"] == tdr.IPR_DOMAIN].drop_duplicates()

    sel = ip[ip["ann"].isin(set(quinasas))]
    return pd.DataFrame([{
        "patron": tdr.PATRON_QUINASA_FUENTE,
        "anotaciones_quinasa": sel["ann"].nunique(),
        "anotaciones_dominio": ip["ann"].nunique(),
        "pct_anotaciones": 100 * sel["ann"].nunique() / ip["ann"].nunique(),
        "proteinas_con_quinasa": sel["target_id"].nunique(),
        "proteinas_con_dominio": ip["target_id"].nunique(),
        "pct_proteinas": 100 * sel["target_id"].nunique() / ip["target_id"].nunique(),
    }])


def bloque_constante(quinasas, por_linea=4):
    """`ANN_QUINASA` lista para pegar en `comun/tdr.py`.

    Se emite ordenada por codigo y con el accession al lado, que es la forma en
    que la constante se puede auditar entrada por entrada.
    """
    items = [f'{c}: "{a}",' for c, a in sorted(quinasas.items())]
    return "\n".join("    " + " ".join(items[i:i + por_linea])
                     for i in range(0, len(items), por_linea))


if __name__ == "__main__":
    entradas = entradas_interpro()
    mapa = _mapa_interpro()
    q = anotaciones_quinasa(entradas, mapa)
    imp = impacto(q)
    print("ANN_QUINASA = {")
    print(bloque_constante(q))
    print("}")
    print(f"\n# {len(q)} anotaciones de quinasa sobre "
          f"{int(imp['anotaciones_dominio'][0])} dominios "
          f"({float(imp['pct_anotaciones'][0]):.2f} %), en "
          f"{int(imp['proteinas_con_quinasa'][0])} de "
          f"{int(imp['proteinas_con_dominio'][0])} proteinas "
          f"({float(imp['pct_proteinas'][0]):.2f} %)")
