"""
figuras_huerfanas.py — las figuras de `huerfanas/`, desde sus tablas.

La corrida de esta carpeta es la más cara del repositorio (carga la capa química
y prioriza una droga por vez); las figuras no la repiten: leen lo que dejó cada
celda de corrida en `resultados/huerfanas_out/`. Se regeneran con

    python comun/figuras.py huerfanas
"""
import numpy as np

import tdr
from figuras import figura

H = "huerfanas"
PSEUDO = "01_*_pseudohuerfanas.csv"
CORTE = 0.1                          # frank < 0.1 = blanco recuperado

# Los organismos que representan a los grupos de `tdr.META` en `03_f06`: de cada
# grupo, el que mas pseudohuerfanas informativas aporta (168 de las 270). Con los
# 16 el panel es una mancha; la seleccion es por grupo, no por resultado.
REPRESENTATIVOS = [3, 1, 21, 17, 10, 26]    # H. sapiens, E. coli, A. thaliana,
                                            # D. melanogaster, S. cerevisiae, P. falciparum


def _etiqueta(sp):
    return tdr.NOMBRE_CORTO.get(sp, sp)


def _pseudohuerfanas(t, con_grupo=False):
    res = t.varias("01", "*_pseudohuerfanas", columna="archivo")
    if con_grupo:
        k1 = t("01_k1_drugs.csv")[["drug_id", "sp_id", "grupo"]]
        res = res.merge(k1.rename(columns={"sp_id": "especie"}),
                        on=["drug_id", "especie"], how="left")
    return res


def _semilla_por_especie(res, plt, titulo):
    """Clase de semilla por especie, en % (informativa = techo del método)."""
    tab = res.groupby(["especie", "semilla"]).size().unstack(fill_value=0)
    for c in ("nula", "no informativa", "informativa"):
        if c not in tab.columns:
            tab[c] = 0
    tab = tab[["nula", "no informativa", "informativa"]]
    tab = tab.div(tab.sum(axis=1), axis=0).sort_values("informativa")
    fig, ax = plt.subplots(figsize=(7, 3.5), tight_layout=True)
    abajo = np.zeros(len(tab))
    for c in tab.columns:
        ax.bar(range(len(tab)), 100 * tab[c], bottom=abajo, color=tdr.COLOR_SEMILLA[c], label=c)
        abajo += 100 * tab[c].values
    ax.set_xticks(range(len(tab)))
    ax.set_xticklabels([_etiqueta(s) for s in tab.index], rotation=90, fontsize=7)
    ax.set_ylabel("% de pseudohuérfanas")
    ax.set_title(titulo)
    ax.legend(fontsize=8)
    return fig


# --- 01 · pseudohuérfanas ---------------------------------------------------

@figura(H, "01_f01_grupos_por_especie",
        verifica='Sin control directo de `01_k1_drugs`; el conjunto evaluado tiene el mismo tamaño que en el oráculo v4.',
        lee="01_k1_drugs.csv")
def grupos_por_especie(t, plt):
    """Pseudohuérfanas por especie, apiladas por grupo 0-3."""
    k1 = t("01_k1_drugs.csv")
    tab = (k1.groupby(["sp_id", "grupo"]).size().unstack(fill_value=0)
           .reindex(columns=[0, 1, 2, 3], fill_value=0))
    tab = tab.loc[tab.sum(axis=1).sort_values(ascending=False).index]
    fig, ax = plt.subplots(figsize=(7, 3.5), tight_layout=True)
    abajo = np.zeros(len(tab))
    etiquetas = {0: "ninguna", 1: "enlace a bioactividad +", 2: "multiespecie", 3: "ambas"}
    for g in [0, 1, 2, 3]:
        ax.bar(range(len(tab)), tab[g], bottom=abajo, color=tdr.COLOR_GRUPO[g],
               label=f"grupo {g} — {etiquetas[g]}")
        abajo += tab[g].values
    ax.set_xticks(range(len(tab)))
    ax.set_xticklabels([_etiqueta(s) for s in tab.index], rotation=90, fontsize=7)
    ax.set_ylabel("drogas pseudohuérfanas")
    ax.legend(fontsize=7)
    return fig


@figura(H, "01_f02_frank_por_grupo",
        verifica='`test_fuga`: la orfanización quita todas las aristas de la droga. Contra el oráculo v4 cambia por el recorte y por la corrección de la semilla (PROVENANCE §3.16), que el oráculo no tiene.',
        lee=(PSEUDO, "01_k1_drugs.csv"))
def frank_por_grupo(t, plt):
    """Distribución de frank por grupo (Fig 2 del paper 2016)."""
    df = _pseudohuerfanas(t, con_grupo=True)
    fig, axes = plt.subplots(4, 1, figsize=(5.5, 6), sharex=True, sharey=True,
                             tight_layout=True)
    bins = np.arange(0, 0.61, 0.01)
    for g, ax in zip([0, 1, 2, 3], axes):
        d = df[df["grupo"] == g]
        ax.hist(d["frank"].dropna(), bins=bins, color=tdr.COLOR_GRUPO[g])
        ax.axvline(CORTE, color=tdr.INK2, lw=1, ls="--")
        ax.set_ylabel(f"grupo {g}", fontsize=8)
        ax.text(0.98, 0.85, f"{d['frank'].notna().sum()} de {len(d)} con semilla · "
                            f"{(d['frank'] < CORTE).mean():.0%} < {CORTE}",
                transform=ax.transAxes, ha="right", fontsize=7, color=tdr.INK2)
    axes[-1].set_xlabel("frank (posición relativa del blanco verdadero)")
    return fig


@figura(H, "01_f03_semilla_por_especie",
        verifica='Consistencia interna: con semilla no informativa el puntaje del blanco es 0 por construcción, y se observa 0 % recuperado. Contra el oráculo v4 cambia por el recorte y por la corrección de la semilla (PROVENANCE §3.16), que el oráculo no tiene.',
        lee=PSEUDO)
def semilla_por_especie(t, plt):
    """Clase de semilla por especie: el techo teórico del método."""
    return _semilla_por_especie(_pseudohuerfanas(t), plt,
                                "Cobertura de semilla: informativa = techo del método")


@figura(H, "01_f04_recuperacion_por_semilla",
        verifica='La misma consistencia interna que `01_f03`, sobre `01_resumen_frank` (comparada contra el oráculo v4).',
        lee="01_resumen_frank.csv")
def recuperacion_por_semilla(t, plt):
    """% de pseudohuérfanas recuperadas (frank < 0.1) por clase de semilla."""
    resumen = t("01_resumen_frank.csv")
    r = resumen[resumen["particion"] == "semilla"].sort_values("pct_recuperadas")
    fig, ax = plt.subplots(figsize=(6, 3.5), tight_layout=True)
    ax.barh(range(len(r)), r["pct_recuperadas"],
            color=[tdr.COLOR_SEMILLA.get(v, tdr.MUTED) for v in r["valor"]])
    ax.set_yticks(range(len(r)))
    ax.set_yticklabels([f"{v} (n = {n:,})".replace(",", " ") for v, n in zip(r["valor"], r["n"])],
                       fontsize=8)
    ax.set_xlabel(f"% con frank < {CORTE}")
    ax.set_title("La recuperación la decide la clase de semilla")
    return fig


# --- 02 · cobertura de semilla ---------------------------------------------

@figura(H, "02_f02_palancas",
        verifica='Sin control. Dos de las tres palancas no pueden rescatar nada en esta base por construcción, así que el 0 es esperado (ver informe).',
        lee="02_palancas.csv")
def palancas(t, plt):
    """Pseudohuérfanas sin semilla informativa que cada palanca rescata."""
    pal = t("02_palancas.csv")
    tab = pal.groupby("palanca")["rescatada"].agg(["sum", "count", "mean"])
    fig, ax = plt.subplots(figsize=(5.5, 3.5), tight_layout=True)
    ax.bar(range(len(tab)), 100 * tab["mean"], color=tdr.CATEGORICOS[:len(tab)])
    ax.set_xticks(range(len(tab)))
    ax.set_xticklabels(tab.index, fontsize=8)
    ax.set_ylabel("% de pseudohuérfanas rescatadas")
    for i, (n, c) in enumerate(zip(tab["sum"], tab["count"])):
        ax.text(i, 100 * n / c, f"{int(n)}/{int(c)}", ha="center", va="bottom", fontsize=7)
    return fig


@figura(H, "02_f03_techo_teorico",
        verifica='Es la fracción informativa de `huerfanas/01` por organismo: coincide con `01_f03` (consistencia interna).',
        lee="02_cobertura_por_especie.csv")
def techo_teorico(t, plt):
    """Fracción con semilla informativa por especie: ninguna reponderación
    puede superarla."""
    cob = t("02_cobertura_por_especie.csv", index_col=0).sort_values("techo")
    fig, ax = plt.subplots(figsize=(6.5, 3.5), tight_layout=True)
    ax.barh(range(len(cob)), 100 * cob["techo"], color=tdr.S3)
    ax.set_yticks(range(len(cob)))
    ax.set_yticklabels([_etiqueta(s) for s in cob.index], fontsize=7)
    ax.set_xlabel("% con semilla informativa (techo del método)")
    return fig


# --- 03 · desorfanización global -------------------------------------------

@figura(H, "03_f01_distribucion_rg",
        verifica='Contra el oráculo v4 (`03_ranking_global`): cambia por el recorte y por la corrección de la semilla (PROVENANCE §3.16), que el oráculo no tiene.',
        lee="03_ranking_global.csv")
def distribucion_rg(t, plt):
    """Posición del blanco verdadero en el ranking global."""
    con_rg = t("03_ranking_global.csv")
    fig, ax = plt.subplots(figsize=(5.5, 3.5), tight_layout=True)
    ax.hist(con_rg["rG"], bins=np.logspace(0, np.log10(con_rg["rG"].max()), 50), color=tdr.S1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("rG (posición en el ranking global)")
    ax.set_ylabel("drogas")
    ax.set_title("Sobre qué se calcula r*G")
    return fig


@figura(H, "03_f02_recuperacion",
        verifica='Contra el oráculo v4: r*G cambia por el recorte y por la corrección de la semilla (PROVENANCE §3.16), que el oráculo no tiene. No hay control independiente del corte; su sensibilidad a k_σ está en `03_rg_sensibilidad.csv`.',
        lee=("03_curva_recuperacion.csv", "03_rg_estrella.csv"))
def recuperacion(t, plt):
    """Fig 3A: ρ(rG) y λ(rG), con el umbral 3σ que define r*G."""
    curva = t("03_curva_recuperacion.csv")
    rg_star = t("03_rg_estrella.csv")
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


@figura(H, "03_f03_rss",
        verifica='Contra el oráculo v4 (`03_rss_distribucion`): cambia por el recorte y por la corrección de la semilla (PROVENANCE §3.16), que el oráculo no tiene.',
        lee="03_rss_distribucion.csv")
def rss(t, plt):
    """Posición del blanco verdadero dentro de su especie, con el acumulado."""
    d = t("03_rss_distribucion.csv")
    fig, ax = plt.subplots(figsize=(5.5, 3.5), tight_layout=True)
    ax.bar(d["posicion"], d["n"], color=tdr.S1, width=1.0)
    ax2 = ax.twinx()
    ax2.plot(d["posicion"], 100 * d["frac_acumulada"], color=tdr.S2)
    ax2.set_ylabel("% acumulado", color=tdr.S2)
    ax2.grid(False)
    ax.set_xscale("log")
    ax.set_xlabel("rSS (posición dentro de la especie)")
    ax.set_ylabel("drogas")
    return fig


@figura(H, "03_f04_directa_indirecta",
        verifica='Contra el oráculo v4 (`03_clases_inferencia`). La partición 68/32 del trabajo de 2016 no se verificó contra el artículo.',
        lee="03_ranking_global.csv")
def directa_indirecta(t, plt):
    """Fig 4b: rG y rSS por clase de inferencia (directa / indirecta)."""
    df = t("03_ranking_global.csv")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.5), tight_layout=True)
    clases = [c for c in ("directa", "indirecta") if c in set(df["clase"])]
    for ax, col, etiqueta in [(ax1, "rG", "rG (ranking global)"),
                              (ax2, "rSS", "rSS (dentro de la especie)")]:
        ax.boxplot([df.loc[df["clase"] == c, col].dropna() for c in clases], showfliers=False)
        ax.set_xticks(range(1, len(clases) + 1))
        ax.set_xticklabels(clases, fontsize=8)
        ax.set_yscale("log")
        ax.set_ylabel(etiqueta)
    n = df["clase"].value_counts(normalize=True)
    ax1.set_title(" · ".join(f"{c}: {100 * n.get(c, 0):.0f} %" for c in clases), fontsize=9)
    return fig


# --- 04 · aplicación a la especie foco -------------------------------------

@figura(H, "04_f01_embudo",
        verifica='Contra el oráculo v4 (`04_sp26_embudo`): menos tratables, explicado por el recorte.',
        lee="04_embudo.csv")
def embudo(t, plt):
    """Compuestos que sobreviven a cada paso del embudo de la especie foco."""
    e = t("04_embudo.csv")
    fig, ax = plt.subplots(figsize=(6, 3), tight_layout=True)
    ax.barh(range(len(e)), e["n"], color=tdr.ORD[:len(e)])
    ax.set_yticks(range(len(e)))
    ax.set_yticklabels(e["paso"], fontsize=8)
    ax.invert_yaxis()
    ax.set_xscale("log")
    for i, n in enumerate(e["n"]):
        ax.text(n * 1.05, i, f"{n:,}", va="center", fontsize=8, color=tdr.INK2)
    ax.set_xlabel("compuestos")
    return fig


@figura(H, "04_f02_rg_sugerencias",
        verifica='Sin sugerencias bajo r*G: no se dibuja.',
        lee=("04_sugerencias.csv", "03_rg_estrella.csv"))
def rg_sugerencias(t, plt):
    """Sin sugerencias bajo r*G no hay nada que dibujar."""
    sug = t("04_sugerencias.csv")
    if sug.empty:
        return None
    rg_star = float(t("03_rg_estrella.csv")["r_g_estrella"].iloc[0])
    fig, ax = plt.subplots(figsize=(5.5, 3.5), tight_layout=True)
    ax.hist(sug["rG"], bins=np.arange(1, rg_star + 2), color=tdr.S1)
    ax.set_xlabel("rG de la sugerencia")
    ax.set_ylabel("pares compuesto-blanco")
    ax.set_title(f"Sugerencias por debajo de r*G = {rg_star:.0f}")
    return fig


@figura(H, "04_f03_familias",
        verifica='Sin sugerencias bajo r*G: no se dibuja.',
        lee="04_familias.csv")
def familias(t, plt):
    """Sin sugerencias no hay familias que contar."""
    fam = t("04_familias.csv")
    if fam.empty:
        return None
    top = fam.sort_values("compuestos", ascending=False).head(15)
    fig, ax = plt.subplots(figsize=(7, 4), tight_layout=True)
    ax.barh(range(len(top)), top["compuestos"], color=tdr.S1)
    ax.set_yticks(range(len(top)))
    # las anotaciones van por código: nombrarlas exige `mapa_interpro`, privado
    ax.set_yticklabels([str(a)[:60] for a in top["anotaciones"]], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("compuestos huérfanos que la reciben")
    return fig


@figura(H, "03_f05_semilla_informativa",
        verifica='Sin control independiente. Referencias internas: el azar en el panel (a); las dos clases se cuentan sobre la misma tabla que `03_f06`.',
        lee=("03_informativa.csv", "03_informativa_topk.csv"))
def semilla_informativa(t, plt):
    """Las pseudohuérfanas con semilla informativa: dónde cae el blanco, más
    allá de `frank < 0.1`."""
    inf = t("03_informativa.csv")
    topk = t("03_informativa_topk.csv")
    color = {"directa": tdr.S1, "indirecta": tdr.S2}
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 4.6), tight_layout=True)

    # (a) recuperacion top-k contra el azar
    a1.plot(topk["k"], 100 * topk["todas"], color=tdr.INK, label=f"todas (n = {len(inf)})")
    for c in ("directa", "indirecta"):
        n = int((inf["clase"] == c).sum())
        a1.plot(topk["k"], 100 * topk[c], color=color[c], label=f"{c} (n = {n})")
    a1.plot(topk["k"], 100 * topk["azar"], color=tdr.MUTED, ls="--", lw=1, label="azar")
    corte = float((CORTE * inf["n_targets"]).median())
    a1.axvline(corte, color=tdr.INK2, ls=":", lw=1)
    a1.text(corte * 0.85, 4, f"frank = {CORTE}\n(≈ posición {corte:,.0f})".replace(",", " "),
            ha="right", fontsize=7, color=tdr.INK2)
    a1.set_xscale("log")
    a1.set_ylim(0, 102)
    a1.set_xlabel("k (posición dentro de la especie)")
    a1.set_ylabel("% con el blanco en el top-k")
    a1.set_title("a · Recuperación top-k")
    a1.legend(fontsize=7, loc="upper left")

    # (b) frank en escala log: cuan por debajo de 0.1 queda. El bin es el doble
    # de ancho que el de la primera version: con 35 bins el histograma era un
    # peine que no se leia proyectado.
    bins = np.logspace(np.log10(inf["frank"].min()) - 0.1, 0, 18)
    for c in ("directa", "indirecta"):
        d = inf.loc[inf["clase"] == c, "frank"]
        a2.hist(d, bins=bins, histtype="step", lw=1.8, color=color[c],
                label=f"{c} (n = {len(d)})")
    a2.set_xscale("log")
    a2.axvline(CORTE, color=tdr.INK2, ls=":", lw=1)
    tope = a2.get_ylim()[1]
    a2.text(CORTE * 0.8, tope * 0.95, f"corte {CORTE}", ha="right",
            fontsize=7, color=tdr.INK2)
    # la mediana va en el grafico, no en la leyenda
    for i, c in enumerate(("directa", "indirecta")):
        m = inf.loc[inf["clase"] == c, "frank"].median()
        a2.axvline(m, color=color[c], ls="--", lw=1.2)
        a2.text(m * 1.15, tope * (0.78 - 0.1 * i), f"mediana {m:.0e}",
                fontsize=7, color=color[c])
    a2.set_xlabel("frank (posición relativa del blanco)")
    a2.set_ylabel("drogas")
    a2.set_title("b · Qué tan por debajo del corte")
    a2.legend(fontsize=7, loc="upper left")

    for ax in (a1, a2):
        ax.set_box_aspect(1)
    return fig


@figura(H, "03_f06_semilla_rss",
        verifica='Sin control independiente. El efecto del tamaño de semilla se comprobó dentro de cada organismo, no sólo en el total; el panel (b) muestra seis de los 16 organismos, elegidos por grupo de `tdr.META` y por cantidad de casos, no por el resultado.',
        lee="03_informativa.csv")
def semilla_rss(t, plt):
    """Qué diluye la señal: tamaño de la semilla y organismo."""
    from scipy.stats import spearmanr
    inf = t("03_informativa.csv")
    color = {"directa": tdr.S1, "indirecta": tdr.S2}
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 4.6), tight_layout=True)

    # (a) tamaño de la semilla contra posicion, con los dos ejes en la misma
    # escala: las dos magnitudes son conteos de proteinas y se comparan.
    for c in ("directa", "indirecta"):
        d = inf[inf["clase"] == c]
        a1.scatter(d["n_semilla"], d["rSS"], s=10, alpha=0.6, color=color[c], label=c)
    rho = spearmanr(inf["n_semilla"], inf["rSS"]).correlation
    lim = (0.8, 1.2 * max(inf["n_semilla"].max(), inf["rSS"].max()))
    a1.set_xscale("log")
    a1.set_yscale("log")
    a1.set_xlim(lim)
    a1.set_ylim(lim)
    a1.set_xlabel("proteínas en la semilla")
    a1.set_ylabel("rSS (posición del blanco)")
    a1.set_title(f"a · Semillas grandes diluyen (Spearman ρ = {rho:.2f})")
    a1.legend(fontsize=7)

    # (b) por organismo, sobre log10(rSS): el violin necesita un eje lineal, asi
    # que se transforma el dato y se etiquetan los ticks como potencias de 10.
    sel = [s for s in REPRESENTATIVOS if (inf["especie"] == s).any()]
    orden = (inf[inf["especie"].isin(sel)].groupby("especie")["rSS"]
             .median().sort_values(ascending=False).index.tolist())
    rng = np.random.default_rng(0)
    datos = [np.log10(inf.loc[inf["especie"] == sp, "rSS"].to_numpy()) for sp in orden]
    v = a2.violinplot(datos, positions=range(len(orden)), orientation="horizontal",
                      widths=0.85, showextrema=False, showmedians=False)
    for cuerpo in v["bodies"]:
        cuerpo.set_facecolor(tdr.GRID)
        cuerpo.set_edgecolor(tdr.AXIS)
        cuerpo.set_alpha(1)
    for i, sp in enumerate(orden):
        d = inf[inf["especie"] == sp]
        y = i + rng.uniform(-0.22, 0.22, len(d))
        a2.scatter(np.log10(d["rSS"]), y, s=10, alpha=0.7,
                   color=[color[c] for c in d["clase"]], zorder=3)
        a2.plot(np.log10(d["rSS"].median()), i, "|", color=tdr.INK, ms=14, mew=2, zorder=4)
    a2.axvline(1, color=tdr.MUTED, ls="--", lw=1)
    a2.set_yticks(range(len(orden)))
    a2.set_yticklabels([f"{_etiqueta(sp)}  n={int((inf['especie'] == sp).sum())}"
                        for sp in orden], fontsize=7)
    a2.set_xticks([0, 1, 2, 3, 4])
    a2.set_xticklabels(["1", "10", "100", "1 000", "10 000"])
    a2.set_xlabel("rSS (posición del blanco);  | = mediana")
    n_sel = int(inf["especie"].isin(sel).sum())
    a2.set_title(f"b · Por organismo (uno por grupo; {n_sel} de {len(inf)} casos)")

    for ax in (a1, a2):
        ax.set_box_aspect(1)
    return fig
