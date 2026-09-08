"""Comparacion de los resultados de este repositorio contra `oraculo_v4/`.

El oraculo son las mismas tablas calculadas por estos mismos notebooks antes de
dos cambios: el recorte de la base (componentes conexas del grafo de fingerprint
sin ningun cluster con bioactividad positiva) y la codificacion a enteros.

La comparacion sirve porque los dos cambios no afectan a todo por igual, y el
criterio es uno solo: **si el analisis toca la capa quimica, el recorte lo mueve;
si no la toca, no lo mueve**. El recorte elimina compuestos, o sea nodos y
aristas de esa capa; no toca proteinas ni anotaciones.

  * `genome_prioritization/` carga `cargar_db(anotaciones=True)`, sin quimica:
    trabaja sobre anotaciones y blancos conocidos -> ESTABLE
  * `analiceDB/` describe la base entera, componentes conexas incluidas -> CAMBIA
  * `huerfanas/` construye la semilla desde el vecindario quimico de cada
    droga (`quimica=True`), y su notebook 03 lee las salidas del 01 -> CAMBIA

Lo estable tiene que salir identico; lo que cambia tiene que cambiar de una
forma que se pueda explicar. Eso es lo que sostiene el resultado ante alguien
que no vio producirlo.

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

# Que se espera de cada analisis, segun toque o no la capa quimica.
ESPERADO = {
    "analiceDB": "cambia",               # describe la base recortada
    "genome_prioritization": "estable",  # anotaciones + druggables, sin quimica
    "huerfanas": "cambia",               # la semilla sale del vecindario quimico
    "verificacion": "estable",
}


def _numericas(df):
    return [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]


def _norm_codigo(serie):
    """`sp26` -> `26`. El oraculo se anonimizo con el prefijo `sp` y las salidas
    nuevas escriben el codigo pelado: sin esto las claves no alinean y toda la
    comparacion sale vacia (que es como se veia antes de arreglarlo)."""
    if serie.dtype == object:
        return serie.astype(str).str.replace(r"^sp(?=-?\d+$)", "", regex=True)
    return serie.astype(str)


def _alinear(nueva, vieja, clave):
    """Alinea dos tablas por `clave`. Devuelve (comunes_n, comunes_v, solo_n, solo_v)."""
    if clave is None or clave not in nueva.columns or clave not in vieja.columns:
        n = min(len(nueva), len(vieja))
        return nueva.head(n).reset_index(drop=True), vieja.head(n).reset_index(drop=True), \
            max(0, len(nueva) - n), max(0, len(vieja) - n)
    n = nueva.copy(); v = vieja.copy()
    n[clave] = _norm_codigo(n[clave])
    v[clave] = _norm_codigo(v[clave])
    n = n.set_index(clave).sort_index()
    v = v.set_index(clave).sort_index()
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


# Fraccion minima de filas en comun para que comparar fila a fila tenga sentido.
MIN_SOLAPE = 0.5


def _resumen_dist(a, b):
    """Compara dos columnas como distribuciones, no fila a fila.

    Hace falta cuando las dos tablas describen el mismo fenomeno sobre conjuntos
    distintos: las muestras de pseudohuerfanas, por ejemplo, se sortean del
    universo k1, que el recorte cambio, asi que apenas comparten filas. Comparar
    por clave ahi da cero y no dice nada; lo que se puede comparar es la forma.
    """
    a = pd.to_numeric(a, errors="coerce").dropna()
    b = pd.to_numeric(b, errors="coerce").dropna()
    if len(a) == 0 or len(b) == 0:
        return {"n": 0, "iguales": float("nan"), "dmax": float("nan"),
                "dmed": float("nan"), "modo": "distribucion"}
    escala = max(abs(a.median()), abs(b.median()), 1e-12)
    d_med = abs(a.median() - b.median())
    d_mean = abs(a.mean() - b.mean())
    # "iguales" pasa a ser cuan cerca estan las medianas, en escala relativa
    cercania = max(0.0, 1.0 - d_med / escala)
    return {"n": int(min(len(a), len(b))), "iguales": float(cercania),
            "dmax": float(max(d_med, d_mean)), "dmed": float(d_med),
            "modo": "distribucion"}


def comparar_tabla(nueva: pd.DataFrame, vieja: pd.DataFrame, clave=None) -> pd.DataFrame:
    """Compara columna por columna dos versiones de la misma tabla.

    Si las dos tablas comparten pocas filas por `clave`, cambia a comparar
    distribuciones y lo declara en la columna `modo`.
    """
    cn, cv, solo_n, solo_v = _alinear(nueva, vieja, clave)
    solape = len(cn) / max(1, min(len(nueva), len(vieja)))
    por_fila = solape >= MIN_SOLAPE
    filas = []
    if por_fila:
        for c in sorted(set(_numericas(cn)) & set(_numericas(cv))):
            r = _resumen_col(cn[c], cv[c]); r.update({"columna": c, "modo": "fila"})
            filas.append(r)
    else:
        for c in sorted(set(_numericas(nueva)) & set(_numericas(vieja))):
            r = _resumen_dist(nueva[c], vieja[c]); r.update({"columna": c})
            filas.append(r)
    out = pd.DataFrame(filas)
    out.attrs.update({"filas_nueva": len(nueva), "filas_vieja": len(vieja),
                      "solo_nueva": solo_n, "solo_vieja": solo_v,
                      "solape": solape, "modo": "fila" if por_fila else "distribucion"})
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
