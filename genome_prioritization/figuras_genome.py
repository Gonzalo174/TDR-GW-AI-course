"""
figuras_genome.py — las figuras de `genome_prioritization/`, desde sus tablas.

Ninguna vuelve a correr el modelo: el barrido (`01_<sp>.csv`) y el ranking de la
especie foco (`02_ranking_<sp>.csv`) los escriben las celdas de corrida. Lo que
se hace acá es lectura y agregación liviana (elegir el óptimo de cada especie
entre 120 filas, por ejemplo). Se regeneran con

    python comun/figuras.py genome_prioritization
"""
import numpy as np
import pandas as pd

import tdr
from figuras import figura

G = "genome_prioritization"
BARRIDO = "01_[0-9]*.csv"          # un CSV por especie; no incluye 01_resumen_*, 01_control_*


def _barrido(t):
    return t.varias("01", "[0-9]*", columna="especie")


def _optimos(barrido):
    import funciones_genome as fg
    return fg.optimos_por_especie(barrido)


def _etiqueta(sp):
    return tdr.NOMBRE_CORTO.get(sp, sp)


# --- 01 · barrido -----------------------------------------------------------

@figura(G, "01_f01_positivos_por_especie",
        verifica='Deriva de la base; el total de druggables lo fija `test_resultados`.',
        lee="01_resumen_especies.csv")
def positivos_por_especie(t, plt):
    """Blancos druggables por especie: los positivos que el barrido recupera."""
    e = t("01_resumen_especies.csv").sort_values("N_druggable")
    fig, ax = plt.subplots(figsize=(7, 3.5), tight_layout=True)
    ax.barh(range(len(e)), e["N_druggable"], color=tdr.S1)
    ax.set_yticks(range(len(e)))
    ax.set_yticklabels([_etiqueta(s) for s in e["sp_id"]], fontsize=8)
    ax.set_xscale("log")
    ax.set_xlabel("blancos druggables (positivos a recuperar)")
    ax.set_title("Positivos por especie")
    return fig


@figura(G, "01_f02_control_v5",
        verifica='Es el control externo (CONVENCIONES §3): `test_los_optimos_reproducen_el_control` exige los 16 óptimos iguales y ΔAUC01 < 1e-9.',
        lee="01_control_v5.csv")
def control_v5(t, plt):
    """AUC01 de esta corrida contra la del control independiente."""
    c = t("01_control_v5.csv")
    fig, ax = plt.subplots(figsize=(4.5, 4.5), tight_layout=True)
    ax.scatter(c["auc01_v5"], c["auc01_v4"], color=tdr.S1, zorder=3)
    lims = [0.5, 1.0]
    ax.plot(lims, lims, color=tdr.MUTED, ls="--", lw=1)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("AUC01 — control independiente")
    ax.set_ylabel("AUC01 — esta corrida")
    ax.set_title("La corrida sobre la base codificada reproduce el control")
    return fig


@figura(G, "01_f03_auc01_por_especie",
        verifica='Barrido idéntico al oráculo v4 (149 de 149 columnas); óptimos iguales al control (`test_resultados`); `test_fuga` (sin fuga en el leave-one-species-out) y `test_metricas` (pAUC y McClish).',
        lee=BARRIDO)
def auc01_por_especie(t, plt):
    """AUC01 en el óptimo de cada especie, bajo leave-one-species-out."""
    o = _optimos(_barrido(t)).sort_values("AUC01")
    fig, ax = plt.subplots(figsize=(7, 3.5), tight_layout=True)
    ax.barh(range(len(o)), o["AUC01"], color=tdr.S1)
    ax.axvline(0.5, color=tdr.MUTED, ls="--", lw=1)
    ax.set_yticks(range(len(o)))
    ax.set_yticklabels([_etiqueta(s) for s in o["especie"]], fontsize=8)
    ax.set_xlim(0.5, 1.0)
    ax.set_xlabel("AUC01 en el óptimo  (0.5 = azar)")
    ax.set_title("Priorización por especie, LOSO")
    return fig


# --- 02 · óptimo y consistencia --------------------------------------------

@figura(G, f"02_f01_roc_{tdr.SP_FOCO}",
        verifica='La AUC01 de la leyenda se recalcula del ranking guardado y coincide con la del óptimo en `02_optimos_por_especie`; `test_metricas`, `test_fuga`.',
        lee=(f"02_ranking_{tdr.SP_FOCO}.csv", "02_optimos_por_especie.csv"))
def roc_foco(t, plt):
    """ROC global y distribución de scores de la especie foco, en su óptimo."""
    from sklearn.metrics import roc_curve
    sp = tdr.SP_FOCO
    rnk = t(f"02_ranking_{sp}.csv")
    opt = t("02_optimos_por_especie.csv")
    a, b, l, g = opt.loc[opt["especie"] == sp, ["alpha", "beta", "lambda_", "gamma"]].values[0]
    l = None if pd.isna(l) else l
    # el ranking de la especie es el global restringido a sus proteinas, con el
    # faltante de score en 0: lo mismo que hace `tdr.ranking_especie`
    sp_rnk = rnk[rnk["de_la_especie"]].assign(score=lambda d: d["score"].fillna(0.0))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), tight_layout=True)
    fpr, tpr, _ = roc_curve(rnk["druggable"], rnk["score"])
    auc_g = tdr.auc_global(rnk["druggable"], rnk["score"])
    auc01 = tdr.auc01_mcclish(sp_rnk["druggable"], sp_rnk["score"])
    ax1.plot(fpr, tpr, color=tdr.S1, label=f"ROC global (AUC = {auc_g:.3f})")
    ax1.plot([0, 1], [0, 1], color=tdr.MUTED, ls="--", lw=1, label="azar")
    ax1.axvline(tdr.FPR_MAX, color=tdr.S2, lw=1, ls=":",
                label=f"FPR = {tdr.FPR_MAX:g} (AUC01 = {auc01:.3f})")
    ax1.set_xlabel("Tasa de falsos positivos")
    ax1.set_ylabel("Tasa de verdaderos positivos")
    ax1.set_title(_etiqueta(sp))
    ax1.legend(loc="lower right", fontsize=8)

    bins = np.logspace(-10, 1, 30)
    for etiqueta, sel, color in [("druggable", 1, tdr.S1), ("no druggable", 0, tdr.S2)]:
        s = rnk.loc[rnk["druggable"] == sel, "score"].dropna()
        ax2.hist(s, bins=bins, density=True, histtype="step", color=color,
                 label=f"{etiqueta} ({len(s)})")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel("NDS")
    ax2.set_ylabel("densidad")
    ax2.set_title(f"alpha={a} beta={b} lambda={l} gamma={g}", fontsize=9)
    ax2.legend(fontsize=8)
    return fig


@figura(G, "02_f02_grilla_beta",
        verifica='El barrido que resume es idéntico al del oráculo v4 (`10_equivalencia`).',
        lee=BARRIDO)
def grilla_beta(t, plt):
    """Efecto de cada parámetro sobre la AUC01; el panel de beta es la
    comparación G'r (beta=0) vs G'rk (beta>0)."""
    barrido = _barrido(t)
    fig, axes = plt.subplots(1, 4, figsize=(11, 3), sharey=True, tight_layout=True)
    for ax, col in zip(axes, ["alpha", "beta", "lambda_", "gamma"]):
        d = barrido.copy()
        # "nan" = modo hibrido, es una categoria mas. Se escribe a mano porque
        # `astype(str)` conserva el faltante en pandas 3 en vez de convertirlo a
        # la cadena "nan", y entonces `sorted` compara str contra float.
        d[col] = d[col].map(lambda v: "nan" if pd.isna(v) else str(v))
        orden = sorted(d[col].unique())
        ax.boxplot([d.loc[d[col] == v, "AUC01"].dropna() for v in orden], showfliers=False)
        ax.set_xticks(range(1, len(orden) + 1))
        ax.set_xticklabels(orden, fontsize=8)
        ax.set_xlabel(col)
        ax.set_title(col, fontsize=9)
    axes[0].set_ylabel("AUC01")
    axes[1].set_title("beta:  0 = G'r,  >0 = G'rk", fontsize=9)
    return fig


@figura(G, "02_f03_plateau",
        verifica='`02_plateau` idéntica a la del oráculo v4.',
        lee="02_plateau.csv")
def plateau(t, plt):
    """Caída relativa de AUC01 al top-K por especie: cuán plano es el óptimo."""
    pl = t("02_plateau.csv").sort_values("caida_top20")
    fig, ax = plt.subplots(figsize=(7, 3.5), tight_layout=True)
    x = np.arange(len(pl))
    for k, color in zip((5, 10, 20), tdr.ORD):
        ax.plot(x, 100 * pl[f"caida_top{k}"], "o-", color=color, ms=3, label=f"top-{k}")
    ax.set_xticks(x)
    ax.set_xticklabels([_etiqueta(s) for s in pl["especie"]], rotation=90, fontsize=7)
    ax.set_ylabel("caida de AUC01 (%)")
    ax.legend(fontsize=8)
    return fig


@figura(G, "02_f04_consistencia_especies",
        verifica='`02_spearman_especies` idéntica a la del oráculo v4.',
        lee="02_spearman_especies.csv")
def consistencia_especies(t, plt):
    """Concordancia (Spearman) entre especies sobre el perfil de la grilla."""
    sp = t("02_spearman_especies.csv", index_col=0)
    sp.columns = [int(c) for c in sp.columns]
    orden = [s for s in tdr.ESPECIES_16 if s in sp.columns]
    m = sp.loc[orden, orden]
    fig, ax = plt.subplots(figsize=(6, 5), tight_layout=True)
    im = ax.imshow(m.values, vmin=-1, vmax=1, cmap="RdBu_r")
    ax.set_xticks(range(len(orden)))
    ax.set_xticklabels([_etiqueta(s) for s in orden], rotation=90, fontsize=7)
    ax.set_yticks(range(len(orden)))
    ax.set_yticklabels([_etiqueta(s) for s in orden], fontsize=7)
    ax.grid(False)
    fig.colorbar(im, ax=ax, label="Spearman sobre AUC01")
    return fig


@figura(G, "02_f05_enriquecimiento_topk",
        verifica='`02_consistencia` idéntica a la del oráculo v4.',
        lee="02_consistencia.csv")
def enriquecimiento_topk(t, plt):
    """Enriquecimiento de cada valor de parámetro en el top-10 de cada especie."""
    frec = t("02_consistencia.csv")
    fig, axes = plt.subplots(1, 4, figsize=(11, 3), sharey=True, tight_layout=True)
    for ax, col in zip(axes, ["alpha", "beta", "lambda_", "gamma"]):
        f = frec[frec["parametro"] == col].sort_values("valor")
        # el faltante de lambda_ es el modo hibrido y se rotula "nan" (ver arriba)
        etiquetas = f["valor"].map(lambda v: "nan" if pd.isna(v) else str(v))
        ax.bar(etiquetas, f["enriquecimiento"], color=tdr.S1)
        ax.axhline(1, color=tdr.MUTED, ls="--", lw=1)
        ax.set_title(col, fontsize=9)
    axes[0].set_ylabel("enriquecimiento en el top-10")
    return fig


# --- 03 · beta y las quinasas ----------------------------------------------

RESUMEN = "03_resumen_beta.csv"


def _orden_resumen(t, col):
    """El resumen ordenado por `col`, con la etiqueta de especie ya resuelta."""
    r = t(RESUMEN).sort_values(col).reset_index(drop=True)
    r["etiqueta"] = [_etiqueta(s) for s in r["especie"]]
    return r


def _dumbbell(ax, r, col_a, col_b, et_a, et_b):
    """Un renglon por especie con los dos valores unidos: hace ver el tamano y el
    signo del cambio de una sola pasada, que es lo que dos barras apiladas
    esconden."""
    y = np.arange(len(r))
    ax.hlines(y, r[col_b], r[col_a], color=tdr.MUTED, lw=1, zorder=1)
    ax.scatter(r[col_b], y, color=tdr.S2, s=28, zorder=3, label=et_b)
    ax.scatter(r[col_a], y, color=tdr.S1, s=28, zorder=3, label=et_a)
    ax.set_yticks(y)
    ax.set_yticklabels(r["etiqueta"], fontsize=7)
    ax.legend(fontsize=8)


@figura(G, "03_f01_auc01_beta",
        verifica='Las AUC01 de beta=1 salen del mismo barrido que `02_optimos_por_especie` y `test_resultados` fija; la de beta=0 se recalcula con la misma metrica (`test_metricas`).',
        lee=RESUMEN)
def auc01_beta(t, plt):
    """AUC01 con la penalizacion por tamano de categoria y sin ella."""
    r = _orden_resumen(t, "auc01_b1")
    fig, ax = plt.subplots(figsize=(7, 4), tight_layout=True)
    _dumbbell(ax, r, "auc01_b1", "auc01_b0", "β = 1  (G′rk)", "β = 0  (G′r)")
    ax.axvline(0.5, color=tdr.MUTED, ls="--", lw=1)
    ax.set_xlabel("AUC01 a (α, λ, γ) fijos  (0.5 = azar)")
    ax.set_title("Apagar la penalización por categoría grande cuesta AUC01")
    return fig


@figura(G, "03_f02_percentil_quinasas",
        verifica='Percentiles del mismo ranking cuya AUC01 verifica `03_f01`; el conjunto de quinasas lo fija `test_resultados` (391 anotaciones, 13 921 proteínas).',
        lee=RESUMEN)
def percentil_quinasas(t, plt):
    """Donde cae la quinasa tipica en el ranking bajo cada beta."""
    r = _orden_resumen(t, "pct_quinasa_mediana_b1")
    fig, ax = plt.subplots(figsize=(7, 4), tight_layout=True)
    _dumbbell(ax, r, "pct_quinasa_mediana_b1", "pct_quinasa_mediana_b0",
              "β = 1", "β = 0")
    ax.set_xlabel("percentil mediano de las quinasas  (0 = cabeza del ranking)")
    ax.set_title("Posición de la quinasa típica")
    return fig


@figura(G, "03_f03_desplazamiento",
        verifica='El contraste es interno: las dos series salen de la misma tabla y del mismo ranking. Si beta moviera el ranking entero por igual, las dos medianas coincidirían.',
        lee="03_desplazamiento.csv")
def desplazamiento(t, plt):
    """Desplazamiento al pasar de beta=0 a beta=1: quinasas contra el resto."""
    d = t("03_desplazamiento.csv")
    orden = (d[d["grupo"] == "quinasa"].sort_values("delta_mediana")["especie"].tolist())
    fig, ax = plt.subplots(figsize=(7.5, 4), tight_layout=True)
    x = np.arange(len(orden))
    for grupo, color, off in [("quinasa", tdr.S1, -0.15), ("resto", tdr.S2, 0.15)]:
        g = d[d["grupo"] == grupo].set_index("especie").loc[orden]
        ax.vlines(x + off, g["delta_q1"], g["delta_q3"], color=color, lw=2, alpha=0.5)
        ax.plot(x + off, g["delta_mediana"], "o", color=color, ms=4, label=grupo)
    ax.axhline(0, color=tdr.MUTED, ls="--", lw=1)
    ax.set_xticks(x)
    ax.set_xticklabels([_etiqueta(s) for s in orden], rotation=90, fontsize=7)
    ax.set_ylabel("Δ percentil  (β=1 − β=0);  > 0 = baja")
    ax.set_title("El desplazamiento que produce β, por grupo (mediana y cuartiles)")
    ax.legend(fontsize=8)
    return fig


@figura(G, "03_f04_enriquecimiento_top100",
        verifica='El esperado es el que daría un ranking indiferente a ser quinasa: `top_k × n_quinasas / n_targets`, con las tres columnas en la misma tabla.',
        lee=RESUMEN)
def enriquecimiento_top100(t, plt):
    """Quinasas en las primeras 100 posiciones, contra lo que daria el azar."""
    r = _orden_resumen(t, "enriq_top100_b1")
    fig, ax = plt.subplots(figsize=(7.5, 4), tight_layout=True)
    x = np.arange(len(r))
    ax.bar(x - 0.2, r["enriq_top100_b0"], width=0.4, color=tdr.S2, label="β = 0")
    ax.bar(x + 0.2, r["enriq_top100_b1"], width=0.4, color=tdr.S1, label="β = 1")
    ax.axhline(1, color=tdr.MUTED, ls="--", lw=1)
    ax.set_xticks(x)
    ax.set_xticklabels(r["etiqueta"], rotation=90, fontsize=7)
    ax.set_ylabel("enriquecimiento de quinasas en el top-100")
    ax.set_title("La cabeza del ranking  (1 = lo que daría el azar)")
    ax.legend(fontsize=8)
    return fig
