"""Comparacion de los resultados de este repositorio contra `oraculo_v4/`.

El oraculo son las mismas tablas calculadas por estos mismos notebooks antes de
dos cambios: el recorte de la base (componentes conexas del grafo de fingerprint
sin ningun cluster con bioactividad positiva) y la codificacion a enteros.

La comparacion sirve porque los dos cambios no afectan a todo por igual:

  * `analiceDB/` describe la base y sus componentes -> el recorte DEBE moverlo
  * `genome_prioritization/` y `huerfanas/` trabajan sobre la capa de
    anotaciones y los druggables -> NO deberian moverse

Lo que debe cambiar cambia de forma explicable y lo que no debe cambiar no
cambia: eso es lo que sostiene el resultado ante alguien que no vio producirlo.

PRINCIPALES
  comparar_tabla(nueva, vieja, clave)   una tabla contra su par del oraculo
  comparar_todo(salidas, oraculo)       recorre las tablas homologas
  clasificar(fila, esperado)            idéntico / explicado / DISCREPANCIA

AUXILIARES
  _alinear, _numericas, _resumen_col
"""
from __future__ import annotations

import pathlib

import numpy as np
import pandas as pd

# Tolerancia relativa para considerar dos numeros "el mismo numero". 1e-9 es
# ruido de punto flotante; por encima de eso hay una diferencia real de calculo.
RTOL = 1e-9

# Que se espera de cada analisis, segun de que depende (ver el docstring).
ESPERADO = {
    "analiceDB": "cambia",              # describe la base recortada
    "genome_prioritization": "estable",  # anotaciones + druggables
    "huerfanas": "estable",
}


def _numericas(df):
    return [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]


def _alinear(nueva, vieja, clave):
    """Alinea dos tablas por `clave`. Devuelve (comunes_n, comunes_v, solo_n, solo_v)."""
    if clave is None or clave not in nueva.columns or clave not in vieja.columns:
        n = min(len(nueva), len(vieja))
        return nueva.head(n).reset_index(drop=True), vieja.head(n).reset_index(drop=True), \
            max(0, len(nueva) - n), max(0, len(vieja) - n)
    n = nueva.set_index(clave).sort_index()
    v = vieja.set_index(clave).sort_index()
    com = n.index.intersection(v.index)
    return (n.loc[com].reset_index(), v.loc[com].reset_index(),
            len(n.index.difference(v.index)), len(v.index.difference(n.index)))


def _resumen_col(a, b):
    """Diferencia entre dos series numericas alineadas."""
    a = pd.to_numeric(a, errors="coerce").to_numpy(dtype=float)
    b = pd.to_numeric(b, errors="coerce").to_numpy(dtype=float)
    ok = ~(np.isnan(a) | np.isnan(b))
    if ok.sum() == 0:
        return {"n": 0, "iguales": np.nan, "dmax": np.nan, "dmed": np.nan}
    a, b = a[ok], b[ok]
    d = np.abs(a - b)
    escala = np.maximum(np.abs(a), np.abs(b))
    rel = np.divide(d, escala, out=np.zeros_like(d), where=escala > 0)
    return {"n": int(ok.sum()),
            "iguales": float((rel <= RTOL).mean()),
            "dmax": float(d.max()),
            "dmed": float(np.median(d))}


def comparar_tabla(nueva: pd.DataFrame, vieja: pd.DataFrame, clave=None) -> pd.DataFrame:
    """Compara columna por columna dos versiones de la misma tabla."""
    cn, cv, solo_n, solo_v = _alinear(nueva, vieja, clave)
    filas = []
    for c in sorted(set(_numericas(cn)) & set(_numericas(cv))):
        r = _resumen_col(cn[c], cv[c])
        r.update({"columna": c})
        filas.append(r)
    out = pd.DataFrame(filas)
    out.attrs.update({"filas_nueva": len(nueva), "filas_vieja": len(vieja),
                      "solo_nueva": solo_n, "solo_vieja": solo_v})
    return out


def clasificar(pct_iguales, analisis) -> str:
    """idéntico / explicado por el recorte / DISCREPANCIA."""
    if pd.isna(pct_iguales):
        return "sin datos"
    if pct_iguales >= 0.999:
        return "idéntico"
    if ESPERADO.get(analisis) == "cambia":
        return "explicado por el recorte"
    return "DISCREPANCIA"


CLAVES = {   # por que columna se alinean las tablas homologas
    "01_disponibilidad_especie.csv": "sp_id",
    "01_dimensiones_capas.csv": None,
    "02_optimos_por_especie.csv": "especie",
    "02_plateau.csv": "especie",
    "02_spearman_especies.csv": "especie",
    "03_bioactividades_por_tag.csv": None,
    "01_resumen_frank.csv": None,
}


def comparar_todo(salidas: pathlib.Path, oraculo: pathlib.Path) -> pd.DataFrame:
    """Recorre las tablas que existen en los dos lados y las compara."""
    filas = []
    for sub in sorted(oraculo.glob("*_out")):
        analisis = sub.name.replace("_out", "")
        for viejo in sorted(sub.glob("*.csv")):
            nuevo = salidas / sub.name / viejo.name
            if not nuevo.exists():
                filas.append({"analisis": analisis, "tabla": viejo.name,
                              "columna": "-", "estado": "sin par nuevo"})
                continue
            try:
                dn, dv = pd.read_csv(nuevo), pd.read_csv(viejo)
            except Exception as e:
                filas.append({"analisis": analisis, "tabla": viejo.name,
                              "columna": "-", "estado": f"ilegible: {type(e).__name__}"})
                continue
            cmp = comparar_tabla(dn, dv, CLAVES.get(viejo.name, "especie"))
            for _, r in cmp.iterrows():
                filas.append({"analisis": analisis, "tabla": viejo.name,
                              "columna": r["columna"], "n": r["n"],
                              "pct_iguales": r["iguales"], "dmax": r["dmax"],
                              "estado": clasificar(r["iguales"], analisis),
                              "filas_nueva": cmp.attrs["filas_nueva"],
                              "filas_vieja": cmp.attrs["filas_vieja"]})
            if cmp.empty:
                filas.append({"analisis": analisis, "tabla": viejo.name,
                              "columna": "(sin columnas numericas comunes)",
                              "estado": "sin datos",
                              "filas_nueva": cmp.attrs["filas_nueva"],
                              "filas_vieja": cmp.attrs["filas_vieja"]})
    return pd.DataFrame(filas)
