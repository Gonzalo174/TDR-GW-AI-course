"""
funciones_genome.py — piezas del barrido de parametros y de la eleccion del optimo.

Es el `full_genome_v5.ipynb` desarmado en funciones (README §1.2): el ciclo que
las combina no vive aca, vive en la celda de corrida del notebook, que es donde
uno lo va a mirar.

Division (README §5):
  PRINCIPAL = cambia un numero que va al paper (semilla, relevance score,
              propagacion, metrica, ranking) -> lleva docstring y test.
  AUXILIAR  = no lo cambia (I/O, paralelizacion, formateo, figuras, chequeos).

Nada se ejecuta al importar este modulo.

Paralelizacion por fork: `evaluar_combinacion` lee `sta` y el contexto de la
especie de globals del modulo, que el proceso padre llena con `fijar_contexto()`
ANTES de crear el pool. Asi las tablas grandes no se vuelven a serializar por
tarea.
"""

from pathlib import Path

import numpy as np
import pandas as pd

import tdr
import nucleo as nf

# Grilla del barrido: 5 x 4 x 4 x 3 = 240 combinaciones nominales, 120 validas
# (cuando lambda_ esta fijo el modelo no usa gamma; ver `grilla_valida`).
PARAMS = {
    "alpha":   [0, 0.3, 0.6, 0.9, 1.2],
    "beta":    [0, 0.5, 1.0, 1.5],
    "lambda_": [0.0, 0.5, 1.0, None],   # None == hibrido (el NA del codigo R)
    "gamma":   [0.8, 1.0, 1.2],
}

# Contexto de la especie en curso; lo llena el padre, lo leen los hijos por fork.
CTX = {}
STA = None


# ============================ PRINCIPALES ============================

def preparar_especie(datos, sp_code, params=PARAMS):
    """LOSO: semilla global sin `sp_code`, cat_rs, true positives y targets de la especie.

    Todo lo que depende solo de la especie —y no de la combinacion de
    parametros— se calcula una vez aca, en el proceso padre.

    Devuelve un dict con
      sp_code     codigo de la especie retenida
      seed        semilla de druggables con peso uniforme, sin evidencia de sp_code
      cat_rs      p-valores de Fisher por categoria (caros: no dependen de alpha)
      tp          positivos a evaluar: druggables de la especie que NO son semilla
      sp_targets  todas las proteinas de la especie (el universo del ranking)
    """
    seed = tdr.semilla_global(datos.posdt, datos.st, sp_out=sp_code)
    cat_rs = nf.get_annot_druggability_pv(datos.sta, seed)

    tp = datos.st[(datos.st["sp_id"] == sp_code) &
                  (datos.st["target_id"].isin(datos.posdt["target_id"]))]
    tp = tp["target_id"].astype(str).unique()
    # los que ya estan en la semilla no se evaluan: serian fuga (ver test_fuga)
    tp = np.setdiff1d(tp, seed["target_id"].astype(str).values)
    sp_targets = datos.st[datos.st["sp_id"] == sp_code]["target_id"].astype(str).unique()

    return {"sp_code": sp_code, "seed": seed, "cat_rs": cat_rs,
            "tp": tp, "sp_targets": sp_targets}


def relevance_por_alpha(ctx, alphas):
    """RS para cada alpha de la grilla (una vez por especie, no por combinacion).

    El test de Fisher es la parte cara del pipeline y no depende de alpha; alpha
    solo reescala los p-valores ya calculados. Devuelve una lista alineada con
    `alphas`, indexable por el `ia` de la combinacion.
    """
    return [tdr.relevance_scores(alpha=a, pv=ctx["cat_rs"])[["ann", "RS"]]
            for a in alphas]


def evaluar_combinacion(combo, params=PARAMS):
    """Una combinacion (ia, ib, il, ig) -> (pcode, AUC01, AUC, Youden, N_targets).

    Se ejecuta en un proceso hijo y lee `STA` y `CTX` como globals (ver el
    encabezado del modulo). `pcode` son los cuatro indices separados por espacio,
    que `guardar_barrido` vuelve a traducir a los valores legibles.

    AUC01 es la pAUC(FPR<=0.10) normalizada y corregida por McClish, calculada
    sobre el ranking restringido a la especie; AUC y el corte de Youden se
    calculan sobre el ranking global. La metrica es la de `comun/tdr.py`: aca no
    se reimplementa (README §1.1).
    """
    ia, ib, il, ig = combo
    beta_val = params["beta"][ib]
    lam_val = params["lambda_"][il]
    gamma_val = params["gamma"][ig]
    pcode = f"{ia + 1} {ib + 1} {il + 1} {ig + 1}"

    rnk = nf.nds(STA, CTX["seed"], lambda_=lam_val, hybrid_gamma=gamma_val,
                 cat_rs=CTX["rs"][ia], beta=beta_val)
    rnk["target_id"] = rnk["target_id"].astype(str)
    rnk["druggable"] = rnk["target_id"].isin(CTX["tp"]).astype(int)

    sp_rnk = rnk[rnk["target_id"].isin(CTX["sp_targets"])].copy()
    sp_rnk["score"] = sp_rnk["score"].fillna(0.0)
    sp_rnk["druggable"] = sp_rnk["druggable"].fillna(0).astype(int)

    auc01 = tdr.auc01_mcclish(sp_rnk["druggable"], sp_rnk["score"])
    auc_g = tdr.auc_global(rnk["druggable"], rnk["score"])
    yo = tdr.youden(rnk["druggable"], rnk["score"])

    return pcode, auc01, auc_g, yo, len(sp_rnk)


def guardar_barrido(filas, salidas, nb, sp_code, params=PARAMS):
    """Arma el DataFrame de resultados, lo ordena por AUC01 y lo escribe como
    `<salidas>/<nb>_<especie>.csv` (README §1.9).

    Traduce el `pcode` de indices a los valores legibles de alpha/beta/lambda_/gamma
    (`lambda_` vacio = modo hibrido, en el que manda gamma).
    """
    df = pd.DataFrame(filas, columns=["params", "AUC01", "AUC", "YoudenCutOff", "N_targets"])
    idx = df["params"].str.split(" ", expand=True).astype(int)
    for j, col in enumerate(["alpha", "beta", "lambda_", "gamma"]):
        df[col] = idx[j].map(lambda i, c=col: params[c][i - 1])
    df = df.sort_values("AUC01", ascending=False).reset_index(drop=True)

    destino = Path(salidas) / f"{nb}_{sp_code}.csv"
    df.to_csv(destino, index=False)
    print(f"  {sp_code}: {len(df)} combinaciones -> {destino.name} "
          f"(AUC01 max {df['AUC01'].iloc[0]:.4f})", flush=True)
    return df


def optimos_por_especie(barrido, k=1):
    """Fila(s) optima(s) por especie; k>1 devuelve el top-K para el plateau.

    `barrido` es lo que devuelve `cargar_barrido`: los CSV del notebook 01
    concatenados con la columna `especie`.
    """
    return (barrido.sort_values(["especie", "AUC01"], ascending=[True, False])
            .groupby("especie", as_index=False, group_keys=False)
            .head(k)
            .reset_index(drop=True))


def consistencia(barrido, k=10):
    """Concordancia entre especies sobre el perfil de las 120 combinaciones.

    Dos vistas, las dos del README §4 (`02_optimo_y_consistencia`):
      * `frecuencias`: cuantas especies ponen cada valor de parametro en su
        top-K, contra la frecuencia marginal de la grilla (enriquecimiento > 1
        = el valor aparece mas de lo que daria el azar);
      * `spearman`: matriz de correlacion de rangos entre los perfiles de AUC01
        de cada par de especies.
    """
    top = optimos_por_especie(barrido, k=k)
    filas = []
    for col in ["alpha", "beta", "lambda_", "gamma"]:
        marg = barrido.groupby("especie")[col].value_counts(normalize=True, dropna=False)
        marg = marg.groupby(level=1, dropna=False).mean()
        obs = top[col].value_counts(normalize=True, dropna=False)
        for val, f in obs.items():
            esperado = marg.get(val, np.nan)
            filas.append({"parametro": col, "valor": val, "frec_topk": f,
                          "frec_grilla": esperado, "enriquecimiento": f / esperado})
    frecuencias = pd.DataFrame(filas)

    perfil = barrido.pivot_table(index="params", columns="especie", values="AUC01")
    spearman = perfil.corr(method="spearman")
    return frecuencias, spearman


def plateau(barrido, tolerancias=(0.001, 0.005, 0.01)):
    """Cuan plano es el entorno del maximo, por especie.

    Si el optimo esta en una meseta ancha, elegir el argmax exacto es arbitrario
    y el numero que va al paper deberia reportarse con esa meseta. Devuelve la
    caida relativa al top-5/10/20 y el tamano del conjunto casi-optimo bajo cada
    tolerancia.
    """
    filas = []
    for sp, g in barrido.groupby("especie"):
        s = g["AUC01"].sort_values(ascending=False).values
        f = {"especie": sp, "n_combos": len(s), "auc01_max": s[0],
             "auc01_min": s[-1]}
        for k in (5, 10, 20):
            f[f"caida_top{k}"] = (s[0] - s[min(k, len(s)) - 1]) / s[0]
        for t in tolerancias:
            f[f"n_dentro_{t}"] = int((s >= s[0] * (1 - t)).sum())
        filas.append(f)
    return pd.DataFrame(filas).sort_values("auc01_max", ascending=False).reset_index(drop=True)


# ============================= AUXILIARES ============================

def grilla_valida(params=PARAMS):
    """(ia, ib, il, ig) sin las repeticiones de gamma.

    Cuando lambda_ esta fijo (no None) el modelo no usa gamma, asi que solo se
    corre el primer gamma: 120 combinaciones en vez de 240.
    """
    combos = []
    for ia, _ in enumerate(params["alpha"]):
        for ib, _ in enumerate(params["beta"]):
            for il, lam_val in enumerate(params["lambda_"]):
                for ig, _ in enumerate(params["gamma"]):
                    if lam_val is not None and ig != 0:
                        continue
                    combos.append((ia, ib, il, ig))
    return combos


def fijar_contexto(sta, ctx):
    """Publica `sta` y el contexto de la especie como globals del modulo, para
    que los hijos del pool los vean por fork. Se llama en el padre, antes de
    `paralelizar`."""
    global STA, CTX
    STA = sta
    CTX = ctx


def paralelizar(func, items, n_core=tdr.N_CORE_DEFAULT, **kw):
    """Importada de `comun/tdr.py`; se reexporta para que la celda de corrida
    hable de un solo modulo."""
    return tdr.paralelizar(func, items, n_core=n_core, **kw)


def cargar_barrido(salidas, nb="01", spoi=None):
    """Concatena los CSV del barrido agregando la columna `especie`."""
    b = tdr.cargar_csvs(salidas, nb, columna="especie")
    return b[b["especie"].isin(spoi)].copy() if spoi else b


def escribir_meta(salidas, nb, **campos):
    """Reexporta `tdr.escribir_meta` (README §1.4)."""
    return tdr.escribir_meta(salidas, nb, **campos)


def comparar_con_control(optimos, control=None):
    """Control del README §1.3: los optimos nuevos tienen que coincidir con los
    de `gon3/resultados/genoma_completo_v5/`.

    Devuelve una tabla por especie con el optimo de cada corrida y la diferencia
    de AUC01. `genoma_completo_v5` ya usa la metrica corregida, asi que la
    columna `coincide_params` deberia dar True en las 16.
    """
    control = Path(control) if control else tdr.CONTROL_V5
    filas = []
    for _, r in optimos.iterrows():
        p = control / f"{r['especie']}.csv"
        if not p.exists():
            filas.append({"especie": r["especie"], "control": "falta"})
            continue
        c = pd.read_csv(p).sort_values("AUC01", ascending=False).iloc[0]
        filas.append({
            "especie": r["especie"],
            "auc01_v4": r["AUC01"], "auc01_v5": c["AUC01"],
            "delta": r["AUC01"] - c["AUC01"],
            "params_v4": _clave(r), "params_v5": _clave(c),
            "coincide_params": _clave(r) == _clave(c),
        })
    return pd.DataFrame(filas)


def _clave(fila):
    """(alpha, beta, lambda_, gamma) como tupla comparable, con NaN -> None."""
    return tuple(None if pd.isna(fila[c]) else float(fila[c])
                 for c in ["alpha", "beta", "lambda_", "gamma"])


def resumen_especies(datos, spoi, optimos=None):
    """Proteinas y druggables por especie; con `optimos`, el optimo al lado.

    Es la tabla que la celda 3 muestra como vista de los datos de entrada: sobre
    que se va a correr el barrido.
    """
    drg = set(datos.posdt["target_id"].astype(str))
    r = (datos.st[datos.st["sp_id"].isin(spoi)]
         .groupby("sp_id")["target_id"].nunique().reset_index(name="N"))
    r["N_druggable"] = r["sp_id"].map(
        lambda s: len(set(datos.st.loc[datos.st["sp_id"] == s, "target_id"].astype(str)) & drg))
    r = r.merge(tdr.META, left_on="sp_id", right_on="sp", how="left").drop(columns=["sp"])
    if optimos is not None:
        r = r.merge(optimos, left_on="sp_id", right_on="especie", how="left")
    return r.sort_values("N_druggable", ascending=False).reset_index(drop=True)


# --- figuras ---------------------------------------------------------------

def fig_roc(sp_rnk, rnk, sp_code, params_opt, plt):
    """ROC global + distribucion de scores de la especie elegida (celda 6 de 02)."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), tight_layout=True)

    from sklearn.metrics import roc_curve
    fpr, tpr, _ = roc_curve(rnk["druggable"], rnk["score"])
    auc_g = tdr.auc_global(rnk["druggable"], rnk["score"])
    auc01 = tdr.auc01_mcclish(sp_rnk["druggable"], sp_rnk["score"])
    ax1.plot(fpr, tpr, color=tdr.S1, label=f"ROC global (AUC = {auc_g:.3f})")
    ax1.plot([0, 1], [0, 1], color=tdr.MUTED, ls="--", lw=1, label="azar")
    ax1.axvline(tdr.FPR_MAX, color=tdr.S2, lw=1, ls=":",
                label=f"FPR = {tdr.FPR_MAX:g} (AUC01 = {auc01:.3f})")
    ax1.set_xlabel("Tasa de falsos positivos")
    ax1.set_ylabel("Tasa de verdaderos positivos")
    ax1.set_title(f"{tdr.NOMBRE_CORTO.get(sp_code, sp_code)}")
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
    a, b, l, g = params_opt
    ax2.set_title(f"alpha={a} beta={b} lambda={l} gamma={g}", fontsize=9)
    ax2.legend(fontsize=8)
    return fig


def fig_grilla_beta(barrido, plt):
    """Efecto de cada parametro sobre la AUC01, con beta destacado.

    El panel de beta es la comparacion G'r (beta=0) vs G'rk (beta>0) que el
    README §1.5 saca de `red_gr_vs_grk/`: ya esta contenida en la grilla.
    """
    fig, axes = plt.subplots(1, 4, figsize=(11, 3), sharey=True, tight_layout=True)
    for ax, col in zip(axes, ["alpha", "beta", "lambda_", "gamma"]):
        datos = barrido.copy()
        # "nan" = modo hibrido, es una categoria mas. Se escribe a mano porque
        # `astype(str)` conserva el faltante en pandas 3 en vez de convertirlo a
        # la cadena "nan", y entonces `sorted` compara str contra float.
        datos[col] = datos[col].map(lambda v: "nan" if pd.isna(v) else str(v))
        orden = sorted(datos[col].unique())
        ax.boxplot([datos.loc[datos[col] == v, "AUC01"].dropna() for v in orden],
                   showfliers=False)
        ax.set_xticks(range(1, len(orden) + 1))
        ax.set_xticklabels(orden, fontsize=8)
        ax.set_xlabel(col)
        ax.set_title(col, fontsize=9)
    axes[0].set_ylabel("AUC01")
    axes[1].set_title("beta:  0 = G'r,  >0 = G'rk", fontsize=9)
    return fig


def fig_consistencia(spearman, plt):
    """Matriz de concordancia entre especies sobre el perfil de la grilla."""
    fig, ax = plt.subplots(figsize=(6, 5), tight_layout=True)
    orden = [s for s in tdr.ESPECIES_16 if s in spearman.columns]
    m = spearman.loc[orden, orden]
    im = ax.imshow(m.values, vmin=-1, vmax=1, cmap="RdBu_r")
    ax.set_xticks(range(len(orden)))
    ax.set_xticklabels([tdr.NOMBRE_CORTO.get(s, s) for s in orden], rotation=90, fontsize=7)
    ax.set_yticks(range(len(orden)))
    ax.set_yticklabels([tdr.NOMBRE_CORTO.get(s, s) for s in orden], fontsize=7)
    ax.grid(False)
    fig.colorbar(im, ax=ax, label="Spearman sobre AUC01")
    return fig


def fig_plateau(pl, plt):
    """Caida relativa al top-K por especie: cuan identificado esta el optimo."""
    fig, ax = plt.subplots(figsize=(7, 3.5), tight_layout=True)
    pl = pl.sort_values("caida_top20")
    x = np.arange(len(pl))
    for k, color in zip((5, 10, 20), tdr.ORD):
        ax.plot(x, 100 * pl[f"caida_top{k}"], "o-", color=color, ms=3, label=f"top-{k}")
    ax.set_xticks(x)
    ax.set_xticklabels([tdr.NOMBRE_CORTO.get(s, s) for s in pl["especie"]],
                       rotation=90, fontsize=7)
    ax.set_ylabel("caida de AUC01 (%)")
    ax.legend(fontsize=8)
    return fig
