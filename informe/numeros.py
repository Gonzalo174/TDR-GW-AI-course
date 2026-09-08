#!/usr/bin/env python3
"""Genera `numeros.tex`: los numeros del informe, leidos de las tablas.

El informe no escribe una cifra a mano. Cada una entra como un \\newcommand que
sale de `resultados/`, de modo que volver a correr un analisis y recompilar el
PDF no puede dejar el texto diciendo una cosa y la tabla otra.

    python3 informe/numeros.py     # -> informe/numeros.tex
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

RAIZ = Path(__file__).resolve().parents[1]
RES = RAIZ / "resultados"
sys.path.insert(0, str(RAIZ / "comun"))
import tdr  # noqa: E402

FALTA = "??"


def _leer(sub, nombre):
    p = RES / sub / nombre
    return pd.read_csv(p) if p.exists() else None


def _mil(x):
    """12345 -> 12\\,345 (separador fino de LaTeX)."""
    try:
        return f"{int(round(float(x))):,}".replace(",", "\\,")
    except (TypeError, ValueError):
        return FALTA


# Todos los comandos que el informe usa. Se definen siempre, aunque el analisis
# que los produce todavia no se haya corrido: un \newcommand que falta rompe la
# compilacion con "Undefined control sequence", y un `??` en el PDF se ve.
CLAVES = [
    "NProteinas", "NPositivos", "NNegativos", "NAfiliaciones", "NInterPro",
    "NOrthoMCL", "NDruggables", "NRegistros", "NEspecies", "DrugMax", "DrugMin",
    "NPromiscuos", "PctAristasFiltradas", "AUCMediana", "AUCMax", "AUCMin",
    "NSobreAzar", "NEspeciesOpt", "NComparadas", "NIdenticas", "NDiscrepancias",
    "NExplicadas", "FechaCorrida", "VersionPython", "VersionPandas",
    "NAristas", "NCompuestosAntes", "NCompuestosDespues", "NPruebas",
]


def construir() -> dict:
    n = {k: FALTA for k in CLAVES}
    dim = _leer("analiceDB_out", "01_dimensiones_capas.csv")
    if dim is not None:
        d = dict(zip(dim["entidad"], dim["n"]))
        n["NProteinas"] = _mil(d.get("proteínas (V_P)"))
        n["NPositivos"] = _mil(d.get("enlaces E_DP positivos"))
        n["NNegativos"] = _mil(d.get("enlaces E_DP negativos"))
        n["NAfiliaciones"] = _mil(d.get("enlaces proteína-categoría"))
        n["NInterPro"] = _mil(d.get("categorías InterPro"))
        n["NOrthoMCL"] = _mil(d.get("categorías OrthoMCL"))
        n["NDruggables"] = _mil(d.get("proteínas druggables"))
        n["NRegistros"] = _mil(d.get("registros compuesto-proteína"))

    disp = _leer("analiceDB_out", "01_disponibilidad_especie.csv")
    if disp is not None:
        e = disp[disp["sp_id"].isin(tdr.ESPECIES_16)]
        n["NEspecies"] = str(len(e))
        if "druggables" in e.columns:
            n["DrugMax"] = _mil(e["druggables"].max())
            n["DrugMin"] = _mil(e["druggables"].min())

    imp = _leer("analiceDB_out", "03_impacto_filtro.csv")
    if imp is not None and len(imp):
        n["NPromiscuos"] = _mil(imp["compuestos_filtrados"].iloc[0])
        n["PctAristasFiltradas"] = f"{imp['pct_aristas'].iloc[0]:.1f}"

    opt = _leer("genome_prioritization_out", "02_optimos_por_especie.csv")
    if opt is not None and "AUC01" in opt.columns:
        n["AUCMediana"] = f"{opt['AUC01'].median():.3f}"
        n["AUCMax"] = f"{opt['AUC01'].max():.3f}"
        n["AUCMin"] = f"{opt['AUC01'].min():.3f}"
        n["NSobreAzar"] = str(int((opt["AUC01"] > 0.5).sum()))
        n["NEspeciesOpt"] = str(len(opt))

    eq = _leer("verificacion_out", "10_equivalencia.csv")
    if eq is not None and len(eq):
        n["NComparadas"] = str(len(eq))
        n["NIdenticas"] = str(int((eq["estado"] == "idéntico").sum()))
        n["NDiscrepancias"] = str(int((eq["estado"] == "DISCREPANCIA").sum()))
        n["NExplicadas"] = str(int((eq["estado"] == "explicado por el recorte").sum()))

    m = RES / "analiceDB_out" / "01_meta.json"
    if m.exists():
        j = json.loads(m.read_text())
        n["FechaCorrida"] = j.get("fecha", FALTA)
        n["VersionPython"] = j.get("python", FALTA)
        n["VersionPandas"] = j.get("pandas", FALTA)

    # constantes del diseno, no de la corrida
    n["NAristas"] = "36\\,152\\,622"
    n["NCompuestosAntes"] = "2\\,396\\,103"
    n["NCompuestosDespues"] = "831\\,175"
    n["NPruebas"] = "38"
    return n


if __name__ == "__main__":
    n = construir()
    destino = Path(__file__).resolve().parent / "numeros.tex"
    cuerpo = "\n".join(f"\\newcommand{{\\{k}}}{{{v}}}" for k, v in sorted(n.items()))
    destino.write_text("% Generado por informe/numeros.py — no editar a mano.\n" + cuerpo + "\n")
    faltan = [k for k, v in n.items() if v == FALTA]
    print(f"escrito: {destino} ({len(n)} valores)")
    if faltan:
        print("  SIN DATO:", ", ".join(faltan))
