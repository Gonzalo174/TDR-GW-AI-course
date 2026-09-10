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
    "NAristas", "NCompuestosAntes", "NCompuestosDespues", "NPruebas", "NFiguras",
    # acondicionamiento de la base (DB/meta.json)
    "PctCompuestosConservados", "NComponentesTotales", "NComponentesConservadas",
    "NBlancosCod", "NEspeciesCod", "NDiccionarios", "NCodigos", "SemillaCod",
    # priorizacion
    "AUCFoco", "AUCGlobalFoco", "NDrugFoco", "NSinPuntajeFoco", "NProtozoos", "NProtozoosArriba", "NBetaPositivo", "EnrBetaUno", "EnrBetaCero",
    "PctTopBetaPositivo", "CaidaTopVeinteMax", "SpearmanMediana", "PctSpearmanPositivo",
    # huerfanas
    "NKUno", "NPseudo", "PctRecuperadas", "PctSinSemilla", "NInformativa",
    "PctRecInformativa", "NNoInformativa", "PctRecNoInformativa", "PctRecGrupoCero",
    "PctRecGrupoTres", "TechoMin", "TechoMax", "NMuestraPalancas", "NRescatadas",
    "NConRG", "RGEstrella", "PctDirecta", "PctIndirecta", "RGDirecta", "RGIndirecta",
    "RSSDirecta", "NInfDirecta", "NInfIndirecta", "PctInfTopUno", "PctInfTopDiez",
    "PctInfTopCien", "PctTopDiezDirecta", "PctTopDiezIndirecta", "AzarTopDiez",
    "RSSMedInfDirecta", "RSSMedInfIndirecta", "RhoSemillaRSS", "PctEmpates",
    "FrankCorteRank", "NInfFueraCorte", "NOrgRho", "NOrgRhoPos", "EmbudoActivos", "EmbudoHuerfanos", "EmbudoTratables", "NSugerencias",
    # contraste con v4
    "NComparadasGP", "NIdenticasGP", "NComparadasHU", "NCambianHU", "NComparadasAN",
    "NCambianAN", "PromiscuosVcuatro", "PctAristasVcuatro", "ControlCoinciden", "ControlDeltaMax",
    "PctRecuperadasVcuatro", "PctSinSemillaVcuatro", "RGEstrellaVcuatro",
]


def _pct(x, dec=1):
    try:
        return f"{float(x):.{dec}f}"
    except (TypeError, ValueError):
        return FALTA


def _leer_oraculo(sub, nombre):
    p = RAIZ / "oraculo_v4" / sub / nombre
    return pd.read_csv(p) if p.exists() else None


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

    _priorizacion(n)
    _huerfanas(n)
    _contraste_v4(n, eq)

    # el acondicionamiento de la base: lo registra DB/meta.json al generarla
    _acondicionamiento(n)
    n["NPruebas"] = _contar_pruebas()
    n["NFiguras"] = _contar_figuras()
    return n


def _acondicionamiento(n):
    p = RAIZ / "DB" / "meta.json"
    if not p.exists():
        return
    m = json.loads(p.read_text())
    r = m["recorte"]
    n["NCompuestosAntes"] = _mil(r["compuestos_totales"])
    n["NCompuestosDespues"] = _mil(r["compuestos_conservados"])
    n["PctCompuestosConservados"] = _pct(100 * r["compuestos_conservados"] / r["compuestos_totales"])
    n["NComponentesTotales"] = _mil(r["componentes_totales"])
    n["NComponentesConservadas"] = _mil(r["componentes_conservadas"])
    aristas = [t["filas"] for t in m["tablas"] if t.get("tabla", "").startswith("01_edges")]
    if aristas:
        n["NAristas"] = _mil(aristas[0])
    d = m["diccionarios"]
    n["NBlancosCod"] = _mil(d["blanco"])
    n["NEspeciesCod"] = _mil(d["especie"])
    n["NDiccionarios"] = str(len(d))
    n["NCodigos"] = _mil(sum(d.values()))
    n["SemillaCod"] = str(m["semilla"])


def _priorizacion(n):
    G = "genome_prioritization_out"
    opt = _leer(G, "02_optimos_por_especie.csv")
    if opt is not None:
        foco = opt[opt["especie"] == tdr.SP_FOCO]
        if len(foco):
            n["AUCFoco"] = f"{foco['AUC01'].iloc[0]:.3f}"
            n["AUCGlobalFoco"] = f"{foco['AUC'].iloc[0]:.3f}"
        n["NBetaPositivo"] = str(int((opt["beta"] > 0).sum()))
        # cuantos de los protozoos ocupan los primeros puestos de AUC01
        m = opt.merge(tdr.META, left_on="especie", right_on="sp", how="left")
        k = int((m["grupo"] == "Protozoos").sum())
        n["NProtozoos"] = str(k)
        n["NProtozoosArriba"] = str(int((m.nlargest(k, "AUC01")["grupo"] == "Protozoos").sum()))
    rk = _leer(G, f"02_ranking_{tdr.SP_FOCO}.csv")
    if rk is not None:
        d = rk[rk["de_la_especie"] & (rk["druggable"] == 1)]
        n["NDrugFoco"] = _mil(len(d))
        n["NSinPuntajeFoco"] = _mil((d["score"].isna() | (d["score"] <= 0)).sum())
    con = _leer(G, "02_consistencia.csv")
    if con is not None:
        b = con[con["parametro"] == "beta"].set_index("valor")
        n["EnrBetaUno"] = f"{b.loc[1.0, 'enriquecimiento']:.2f}"
        n["EnrBetaCero"] = f"{b.loc[0.0, 'enriquecimiento']:.2f}"
        n["PctTopBetaPositivo"] = _pct(100 * (1 - b.loc[0.0, "frec_topk"]), 0)
    pl = _leer(G, "02_plateau.csv")
    if pl is not None:
        n["CaidaTopVeinteMax"] = _pct(100 * pl["caida_top20"].max())
    sp = _leer(G, "02_spearman_especies.csv")
    if sp is not None:
        import numpy as np
        m = sp.set_index(sp.columns[0]).values
        fuera = m[~np.eye(len(m), dtype=bool)]
        n["SpearmanMediana"] = f"{np.median(fuera):.2f}"
        n["PctSpearmanPositivo"] = _pct(100 * (fuera > 0).mean(), 0)


def _huerfanas(n):
    H = "huerfanas_out"
    k1 = _leer(H, "01_k1_drugs.csv")
    if k1 is not None:
        n["NKUno"] = _mil(len(k1))
    r = _leer(H, "01_resumen_frank.csv")
    if r is not None:
        def fila(part, valor):
            f = r[(r["particion"] == part) & (r["valor"].astype(str) == str(valor))]
            return f.iloc[0] if len(f) else None
        tot = fila("total", "todas")
        if tot is not None:
            n["NPseudo"] = _mil(tot["n"])
            n["PctRecuperadas"] = _pct(tot["pct_recuperadas"])
            n["PctSinSemilla"] = _pct(tot["pct_sin_semilla"], 0)
        for clave, valor in (("Informativa", "informativa"), ("NoInformativa", "no informativa")):
            f = fila("semilla", valor)
            if f is not None:
                n[f"N{clave}"] = _mil(f["n"])
                n[f"PctRec{clave}"] = _pct(f["pct_recuperadas"])
        for clave, g in (("Cero", 0), ("Tres", 3)):
            f = fila("grupo", g)
            if f is not None:
                n[f"PctRecGrupo{clave}"] = _pct(f["pct_recuperadas"])
    cob = _leer(H, "02_cobertura_por_especie.csv")
    if cob is not None:
        n["TechoMin"] = _pct(100 * cob["techo"].min())
        n["TechoMax"] = _pct(100 * cob["techo"].max(), 0)
    pal = _leer(H, "02_palancas.csv")
    if pal is not None:
        n["NMuestraPalancas"] = _mil(pal.groupby("palanca").size().max())
        n["NRescatadas"] = _mil(pal["rescatada"].sum())
    rg = _leer(H, "03_rg_estrella.csv")
    if rg is not None:
        n["RGEstrella"] = _mil(rg["r_g_estrella"].iloc[0])
    rk = _leer(H, "03_ranking_global.csv")
    if rk is not None:
        n["NConRG"] = _mil(len(rk))
        frac = rk["clase"].value_counts(normalize=True)
        n["PctDirecta"] = _pct(100 * frac.get("directa", 0), 0)
        n["PctIndirecta"] = _pct(100 * frac.get("indirecta", 0), 0)
        med = rk.groupby("clase")[["rG", "rSS"]].median()
        if "directa" in med.index:
            n["RGDirecta"] = _mil(med.loc["directa", "rG"])
            n["RSSDirecta"] = _mil(med.loc["directa", "rSS"])
        if "indirecta" in med.index:
            n["RGIndirecta"] = _mil(med.loc["indirecta", "rG"])
    inf = _leer(H, "03_informativa.csv")
    topk = _leer(H, "03_informativa_topk.csv")
    if inf is not None and topk is not None:
        from scipy.stats import spearmanr
        tk = topk.set_index("k")
        n["NInfDirecta"] = _mil((inf["clase"] == "directa").sum())
        n["NInfIndirecta"] = _mil((inf["clase"] == "indirecta").sum())
        for clave, k in (("Uno", 1), ("Diez", 10), ("Cien", 100)):
            n[f"PctInfTop{clave}"] = _pct(100 * tk.loc[k, "todas"], 0)
        n["PctTopDiezDirecta"] = _pct(100 * tk.loc[10, "directa"], 0)
        n["PctTopDiezIndirecta"] = _pct(100 * tk.loc[10, "indirecta"], 0)
        n["AzarTopDiez"] = _pct(100 * tk.loc[10, "azar"], 2)
        med = inf.groupby("clase")["rSS"].median()
        n["RSSMedInfDirecta"] = f"{med.get('directa'):g}"
        n["RSSMedInfIndirecta"] = f"{med.get('indirecta'):g}"
        n["RhoSemillaRSS"] = f"{spearmanr(inf['n_semilla'], inf['rSS']).correlation:.2f}"
        n["PctEmpates"] = _pct(100 * inf["empate"].mean(), 0)
        # la dilucion, dentro de cada organismo con casos suficientes
        rhos = [spearmanr(d["n_semilla"], d["rSS"]).correlation
                for _, d in inf.groupby("especie") if len(d) >= 15]
        n["NOrgRho"] = str(len(rhos))
        n["NOrgRhoPos"] = str(sum(r > 0 for r in rhos))
        n["FrankCorteRank"] = _mil((0.1 * inf["n_targets"]).median())
        n["NInfFueraCorte"] = _mil((inf["frank"] >= 0.1).sum())
    emb = _leer(H, "04_embudo.csv")
    if emb is not None and len(emb) == 3:
        n["EmbudoActivos"], n["EmbudoHuerfanos"], n["EmbudoTratables"] = (
            _mil(v) for v in emb["n"])
    p = RES / H / "04_sugerencias.csv"
    if p.exists():
        try:
            n["NSugerencias"] = _mil(len(pd.read_csv(p)))
        except pd.errors.EmptyDataError:
            n["NSugerencias"] = "0"


def _contraste_v4(n, eq):
    if eq is not None and len(eq):
        cambia = "explicado por el recorte"
        for clave, analisis in (("GP", "genome_prioritization"), ("HU", "huerfanas"),
                                ("AN", "analiceDB")):
            e = eq[(eq["analisis"] == analisis) & eq["estado"].isin(["idéntico", cambia])]
            n[f"NComparadas{clave}"] = str(len(e))
            if clave == "GP":
                n["NIdenticasGP"] = str(int((e["estado"] == "idéntico").sum()))
            else:
                n[f"NCambian{clave}"] = str(int((e["estado"] == cambia).sum()))
    imp = _leer_oraculo("analiceDB_out", "03_impacto_filtro.csv")
    if imp is not None:
        n["PromiscuosVcuatro"] = _mil(imp["compuestos_filtrados"].iloc[0])
        n["PctAristasVcuatro"] = _pct(imp["pct_aristas"].iloc[0], 0)
    ctl = _leer("genome_prioritization_out", "01_control_v5.csv")
    if ctl is not None:
        n["ControlCoinciden"] = str(int(ctl["coincide_params"].sum()))
        n["ControlDeltaMax"] = f"{ctl['delta'].abs().max():.0e}".replace("0e+00", "0")
    r = _leer_oraculo("huerfanas_out", "01_resumen_frank.csv")
    if r is not None:
        tot = r[r["particion"] == "total"].iloc[0]
        n["PctRecuperadasVcuatro"] = _pct(tot["pct_recuperadas"])
        n["PctSinSemillaVcuatro"] = _pct(tot["pct_sin_semilla"], 0)
    rg = _leer_oraculo("huerfanas_out", "03_rg_estrella.csv")
    if rg is not None:
        n["RGEstrellaVcuatro"] = _mil(rg["r_g_estrella"].iloc[0])


def _contar_figuras():
    """Figuras registradas en `comun/figuras.py`: las que se regeneran desde
    las tablas sin tocar la base."""
    try:
        sys.path.insert(0, str(RAIZ / "comun"))
        import figuras
        return str(len(figuras.cargar_registro()))
    except Exception:
        return FALTA


def _contar_pruebas():
    """Cuenta las pruebas de `comun/tests/` en vez de tenerlas escritas: un
    numero a mano se desactualiza en cuanto se agrega un test, que es
    exactamente lo que este informe dice que no hace."""
    import unittest
    try:
        suite = unittest.defaultTestLoader.discover(str(RAIZ / "comun" / "tests"),
                                                    pattern="test_*.py",
                                                    top_level_dir=str(RAIZ / "comun" / "tests"))
        return str(suite.countTestCases())
    except Exception:
        return FALTA


def tabla_optimos() -> str:
    """El cuadro de optimos por especie, en LaTeX, leido de la tabla."""
    opt = _leer("genome_prioritization_out", "02_optimos_por_especie.csv")
    if opt is None:
        return "\\textbf{??} falta correr genome\\_prioritization/02\n"
    opt = opt.merge(tdr.META, left_on="especie", right_on="sp", how="left")
    filas = []
    for _, r in opt.sort_values("AUC01", ascending=False).iterrows():
        lam = "híbrido" if pd.isna(r["lambda_"]) else f"{r['lambda_']:g}"
        par = "sí" if r.get("parasito") else "no"
        filas.append(f"sp{int(r['especie']):02d} & {r.get('grupo', '')} & {par} & "
                     f"{r['AUC01']:.3f} & {r['AUC']:.3f} & {r['alpha']:g} & {r['beta']:g} & "
                     f"{lam} & {r['gamma']:g} \\\\")
    return ("% Generado por informe/numeros.py — no editar a mano.\n"
            "\\begin{tabular}{@{}lllrrrrrr@{}}\n\\toprule\n"
            "especie & grupo & parásito & AUC01 & AUC & $\\alpha$ & $\\beta$ & "
            "$\\lambda$ & $\\gamma$ \\\\\n\\midrule\n"
            + "\n".join(filas) + "\n\\bottomrule\n\\end{tabular}\n")


if __name__ == "__main__":
    (Path(__file__).resolve().parent / "tabla_optimos.tex").write_text(tabla_optimos())
    n = construir()
    destino = Path(__file__).resolve().parent / "numeros.tex"
    cuerpo = "\n".join(f"\\newcommand{{\\{k}}}{{{v}}}" for k, v in sorted(n.items()))
    destino.write_text("% Generado por informe/numeros.py — no editar a mano.\n" + cuerpo + "\n")
    faltan = [k for k, v in n.items() if v == FALTA]
    print(f"escrito: {destino} ({len(n)} valores)")
    if faltan:
        print("  SIN DATO:", ", ".join(faltan))
