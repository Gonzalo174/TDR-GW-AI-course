"""
nucleo.py
---------
Nucleo del modelo: `nds_fun.py` de la versión v3 del proyecto (fuera de este repositorio) movido tal cual (traduccion de
nds_fun.R). Congelado: se cambia solo con un test de `comun/tests/` que lo
justifique. Todo lo demas —rutas, carga, metricas, paralelizacion, figuras—
vive en `comun/tdr.py`.

Cambios respecto del original (cada uno con su test o su registro):
  2026-09-08  clave de salida con map(str, sp_out) (PROVENANCE.md §3.4).
  2026-09-24  get_druggable_targets consulta ddt/dds con cluster_id y no con
              drug_id (comun/tests/test_vecinos.py, PROVENANCE.md §3.16). Un
              unico salto, como antes. Portado desde TDR_2026_v4.

Dependencies:
    pandas, numpy, scipy

All tables that were data.table in R are now pandas DataFrames.
Column names are preserved exactly as in the original R code.
"""

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import fisher_exact


# ---------------------------------------------------------------------------
# getDruggableTargets
# ---------------------------------------------------------------------------

def get_druggable_targets(
    posdt,
    st,
    doi=None,
    negdt=None,
    use_neighbors=False,
    dds=None,
    ddt=None,
    tclus=None,
    sclus=None,
    tani_exp=1,
    sub_exp=1,
    sp_out=None,
    sp_out_kfold=1,
    rseed=None,
    b_max_consolidate=True,
):
    """
    Build a resource vector defining druggable targets to be used as seeds in
    prioritization tasks.
 
    Parameters
    ----------
    posdt         : DataFrame  drug-target positive bioactivities  [drug_id, target_id, ...]
    st            : DataFrame  target-species table                [target_id, sp_id, ...]
    doi           : list/array drugs of interest (all if None)
    negdt         : DataFrame  drug-target negative bioactivities
    use_neighbors : bool       if True and doi is not None, expand via chemical space
    dds           : DataFrame  drug-drug substructure relationships  (keyed on 'from'/'to')
    ddt           : DataFrame  drug-drug tanimoto similarities       [clusID1, clusID2, weight]
    tclus         : DataFrame  drug-tanimoto-cluster table           [cluster_id, drug]
    sclus         : DataFrame  drug-substructure-cluster table       [cluster_id, drug]
    tani_exp      : float      exponent for tanimoto similarity weights
    sub_exp       : float      exponent for substructure similarity weights
    sp_out        : str/list   species code(s) whose evidence to remove
    sp_out_kfold  : float      fraction of druggable targets to remove (1 = all)
    rseed         : int        random seed (used when sp_out_kfold < 1)
    b_max_consolidate : bool   if True, report the max weight across evidence types per target
 
    Returns
    -------
    dict  { label : DataFrame([target_id, w]) }
          or { label : None } when no evidence is found
    """
 
    if rseed is not None:
        np.random.seed(rseed)
 
    lres = {}
 
    # Join posdt with st on target_id
    sdt = st.merge(posdt, on="target_id", how="right")
    sdt = sdt[sdt["sp_id"].notna()].copy()
 
    # ---- Remove evidence from selected species ----------------------------
    if sp_out is not None:
        if isinstance(sp_out, str):
            sp_out = [sp_out]
        dtout = sdt.loc[sdt["sp_id"].isin(sp_out), "target_id"].unique()
        if sp_out_kfold < 1:
            n_remove = int(len(dtout) * sp_out_kfold)
            dtout = np.random.choice(dtout, n_remove, replace=False)
        sdt = sdt[~((sdt["target_id"].isin(dtout)) & (sdt["sp_id"].isin(sp_out)))].copy()
 
    # ---- No drugs of interest: aggregate all targets ---------------------
    if doi is None:
        a = sdt[["target_id"]].copy()
        a["w"] = 1
        aux = "all"
        if sp_out is not None:
            # map(str, ...): los codigos de especie son enteros en esta base.
            # Arma la clave del dict de salida y no interviene en ningun
            # calculo (PROVENANCE.md §3.4).
            aux += "_sp_out_" + "_".join(map(str, sp_out))
        if sp_out_kfold < 1:
            aux += f"_sp_out_kfold_{sp_out_kfold}"
        lres[aux] = a.groupby("target_id", as_index=False)["w"].sum()
 
    # ---- Drugs of interest, no neighbors ---------------------------------
    else:
        if not use_neighbors:
            sdt_doi = sdt[sdt["drug_id"].isin(doi)].copy()
            if len(sdt_doi) > 0:
                a = sdt_doi[["target_id"]].copy()
                a["w"] = 1
                lres["doi"] = a.groupby("target_id", as_index=False)["w"].sum()
            else:
                lres["doi"] = None
 
        # ---- Drugs of interest WITH neighbors (chemical space) -----------
        else:
            for idoi in doi:
                doi_posdt = posdt.loc[posdt["drug_id"] == idoi, "target_id"].values
                doi_negdt = negdt.loc[negdt["drug_id"] == idoi, "target_id"].values if negdt is not None else np.array([])
 
                # --- Same substructure cluster (sclus) ---
                ssc_clusters = sclus.loc[sclus["drug"] == idoi, "cluster_id"].values
                ssc = sclus.loc[sclus["cluster_id"].isin(ssc_clusters), "drug"].values
                ssc = ssc[ssc != idoi]
                out1 = posdt.loc[posdt["drug_id"].isin(ssc) & posdt["target_id"].isin(doi_negdt), "drug_id"].values
                out2 = negdt.loc[negdt["drug_id"].isin(ssc) & negdt["target_id"].isin(doi_posdt), "drug_id"].values if negdt is not None else np.array([])
                ssc = np.setdiff1d(ssc, np.union1d(out1, out2))
 
                # --- Same tanimoto cluster (tclus) ---
                stc_clusters = tclus.loc[tclus["drug"] == idoi, "cluster_id"].values
                stc = tclus.loc[tclus["cluster_id"].isin(stc_clusters), "drug"].values
                stc = stc[stc != idoi]
                out1 = posdt.loc[posdt["drug_id"].isin(stc) & posdt["target_id"].isin(doi_negdt), "drug_id"].values
                out2 = negdt.loc[negdt["drug_id"].isin(stc) & negdt["target_id"].isin(doi_posdt), "drug_id"].values if negdt is not None else np.array([])
                stc = np.setdiff1d(stc, np.union1d(out1, out2))
 
                # --- Tanimoto clusters of substructure neighbors (stc2) ---
                stc2_clusters = tclus.loc[tclus["drug"].isin(ssc), "cluster_id"].values
                stc2 = tclus.loc[tclus["cluster_id"].isin(stc2_clusters), "drug"].values
                stc2 = stc2[stc2 != idoi]
                out1 = posdt.loc[posdt["drug_id"].isin(stc2) & posdt["target_id"].isin(doi_negdt), "drug_id"].values
                out2 = negdt.loc[negdt["drug_id"].isin(stc2) & negdt["target_id"].isin(doi_posdt), "drug_id"].values if negdt is not None else np.array([])
                stc2 = np.setdiff1d(stc2, np.union1d(out1, out2))
 
                # --- Substructure clusters of tanimoto neighbors (ssc2) ---
                ssc2_clusters = sclus.loc[sclus["drug"].isin(stc), "cluster_id"].values
                ssc2 = sclus.loc[sclus["cluster_id"].isin(ssc2_clusters), "drug"].values
                ssc2 = ssc2[ssc2 != idoi]
                out1 = posdt.loc[posdt["drug_id"].isin(ssc2) & posdt["target_id"].isin(doi_negdt), "drug_id"].values
                out2 = negdt.loc[negdt["drug_id"].isin(ssc2) & negdt["target_id"].isin(doi_posdt), "drug_id"].values if negdt is not None else np.array([])
                ssc2 = np.setdiff1d(ssc2, np.union1d(out1, out2))
 
                # Identity neighbors: union of all cluster-based neighbors
                nei1 = np.unique(np.concatenate([ssc, ssc2, stc, stc2]))
                target1 = posdt.loc[posdt["drug_id"].isin(nei1), ["target_id"]].copy()
                target1["w"] = 1
                target1["type"] = "identity"
                b_target1 = len(target1) > 0
 
                # --- Tanimoto similarity navigation ---
                doi_and_nei = np.concatenate([nei1, [idoi]]).astype(float)
                # ddt y dds unen CLUSTERS identitarios (createDB/01-02): se consultan
                # con los cluster_id de estas drogas, no con sus drug_id. Hasta el
                # 2026-09-24 se consultaban con drug_id (heredado de nds_fun.R, cuando
                # las aristas eran droga-droga) y devolvian los vecinos de otro
                # cluster. Test: comun/tests/test_vecinos.py
                t_cl = tclus.loc[tclus["drug"].isin(doi_and_nei), "cluster_id"].unique()
                s_cl = sclus.loc[sclus["drug"].isin(doi_and_nei), "cluster_id"].unique()
 
                neitc_1 = ddt[ddt["clusID1"].isin(t_cl)][["clusID2", "weight"]].copy()
                neitc_1 = neitc_1.dropna(subset=["weight"])
                neitc_1 = neitc_1.rename(columns={"clusID2": "cluster_id"})
 
                neitc_2 = ddt[ddt["clusID2"].isin(t_cl)][["clusID1", "weight"]].copy()
                neitc_2 = neitc_2.dropna(subset=["weight"])
                neitc_2 = neitc_2.rename(columns={"clusID1": "cluster_id"})
 
                neitc = pd.concat([neitc_1, neitc_2], ignore_index=True)
                cdw = tclus.merge(neitc, on="cluster_id", how="inner")
 
                out1 = posdt.loc[posdt["drug_id"].isin(cdw["drug"]) & posdt["target_id"].isin(doi_negdt), "drug_id"].values
                out2 = negdt.loc[negdt["drug_id"].isin(cdw["drug"]) & negdt["target_id"].isin(doi_posdt), "drug_id"].values if negdt is not None else np.array([])
                t_nei_drugs = np.setdiff1d(cdw["drug"].values, np.union1d(out1, out2))
 
                cdw = cdw[cdw["drug"].isin(t_nei_drugs)]
                cdwt = cdw.merge(posdt.rename(columns={"drug_id": "drug"}), on="drug", how="inner")
 
                b_target_tani = len(cdwt) > 1
                if b_target_tani:
                    target_tani = (
                        cdwt.groupby(["cluster_id", "target_id"])["weight"]
                        .apply(lambda x: (x ** tani_exp).max())
                        .reset_index(name="w")
                    )
                    target_tani = target_tani.groupby("target_id", as_index=False)["w"].sum()
                    target_tani["type"] = "tani"
 
                # --- Substructure similarity navigation ---
                neisc = dds[dds.index.isin(s_cl)][["to", "weight"]].copy()
                neisc = neisc.dropna(subset=["weight"])
                neisc = neisc.rename(columns={"to": "cluster_id"})
 
                cdw = sclus.merge(neisc, on="cluster_id", how="inner")
 
                out1 = posdt.loc[posdt["drug_id"].isin(cdw["drug"]) & posdt["target_id"].isin(doi_negdt), "drug_id"].values
                out2 = negdt.loc[negdt["drug_id"].isin(cdw["drug"]) & negdt["target_id"].isin(doi_posdt), "drug_id"].values if negdt is not None else np.array([])
                s_nei_drugs = np.setdiff1d(cdw["drug"].values, np.union1d(out1, out2))
 
                cdw = cdw[cdw["drug"].isin(s_nei_drugs)]
                cdwt = cdw.merge(posdt.rename(columns={"drug_id": "drug"}), on="drug", how="inner")
 
                b_target_sub = len(cdwt) > 1
                if b_target_sub:
                    target_sub = (
                        cdwt.groupby(["cluster_id", "target_id"])["weight"]
                        .apply(lambda x: (x ** sub_exp).max())
                        .reset_index(name="w")
                    )
                    target_sub = target_sub.groupby("target_id", as_index=False)["w"].sum()
                    target_sub["type"] = "sub"
 
                # --- Combine all evidence ---
                pieces = []
                if b_target1:
                    pieces.append(target1.copy())
                if b_target_tani:
                    pieces.append(target_tani)
                if b_target_sub:
                    pieces.append(target_sub)
 
                sdt_doi = sdt[sdt["drug_id"] == idoi]
                if len(sdt_doi) > 0:
                    direct = sdt_doi[["target_id"]].copy()
                    direct["w"] = 1
                    direct["type"] = "direct"
                    pieces.append(direct)
 
                if pieces:
                    # Ensure all pieces have the same columns before concat
                    for p in pieces:
                        for col in ["target_id", "w", "type"]:
                            if col not in p.columns:
                                p[col] = pd.NA
                    a = pd.concat(pieces, ignore_index=True)
                else:
                    a = pd.DataFrame(columns=["target_id", "w", "type"])
 
                if len(a) > 1:
                    if b_max_consolidate:
                        lres[str(idoi)] = a.groupby("target_id", as_index=False)["w"].max()
                    else:
                        lres[str(idoi)] = a
                else:
                    lres[str(idoi)] = None
 
    return lres
 
 


# ---------------------------------------------------------------------------
# getAnnotDruggabilityPv
# ---------------------------------------------------------------------------

def get_annot_druggability_pv(sta, druggables):
    """
    Calculate the p-value of over-representation of druggable targets
    annotated to a given category (OG, PF), using a one-sided Fisher exact test.

    Parameters
    ----------
    sta         : DataFrame  species-target-annotation table  [target_id, ann, ...]
    druggables  : DataFrame  druggable targets                [target_id, w]

    Returns
    -------
    DataFrame  [ann, nTargets, nDruggables, pv]
    """

    druggable_ids = druggables["target_id"].values

    # targets annotated to each category
    ann_trg = sta.groupby("ann").size().reset_index(name="count")
    # druggables annotated to each category
    ann_drg = (
        sta[sta["target_id"].isin(druggable_ids)]
        .groupby("ann")
        .size()
        .reset_index(name="count")
    )

    a = ann_trg.merge(ann_drg, on="ann", suffixes=("_targets", "_druggables"), how="inner")
    a = a[a["ann"].notna()].copy()
    a = a.rename(columns={"count_targets": "nTargets", "count_druggables": "nDruggables"})

    # Annotation type indices
#     items = {
#         "IP": sta[sta["ann"].str.startswith("IP", na=False)].index,
#         "OG": sta[sta["ann"].str.startswith("OG", na=False)].index,
#     }
    items = {
        "IP": sta[sta["db"] == "ip" ].index,
        "OG": sta[sta["db"] == "omcl" ].index,
    }

    def do_fisher(ab, anb, tot_druggables, tot_non_druggables):
        m = np.array([
            [ab,       tot_druggables - ab],
            [anb, tot_non_druggables - anb],
        ])
        _, pval = fisher_exact(m, alternative="greater")
        return pval

    a["pv"] = np.nan

    for prefix, idx in items.items():
        tot_druggables = sta[sta["target_id"].isin(druggable_ids)]["target_id"].nunique()
        tot_non_druggables = sta.loc[idx, "target_id"].nunique() - tot_druggables

        mask = a["ann"].str.startswith(prefix, na=False)
        for i, row in a[mask].iterrows():
            pv = do_fisher(
                row["nDruggables"],
                row["nTargets"] - row["nDruggables"],
                tot_druggables,
                tot_non_druggables,
            )
            a.at[i, "pv"] = pv

    return a.sort_values("pv").reset_index(drop=True)


# ---------------------------------------------------------------------------
# RS  (Relevance Score)
# ---------------------------------------------------------------------------

def rs(pv, adjust=False, umbral=0.2, alpha=0.6, method="fdr_bh"):
    """
    Estimate a Relevance Score from a vector of p-values.

    Parameters
    ----------
    pv      : array-like  p-values (may be a named Series or plain array)
    adjust  : bool        if True, work with FDR-adjusted q-values
    umbral  : float       if adjust=False: quantile defining the max p-value with RS=1
                          if adjust=True:  q-value cutoff for RS=1
    alpha   : float       exponent applied to (log(p)/log(quantile))
    method  : str         multiple-testing correction method passed to
                          statsmodels (default "fdr_bh" ≈ R's "fdr")

    Returns
    -------
    numpy array (same length as pv) with relevance scores; preserves index if
    pv is a Series.
    """
    from statsmodels.stats.multitest import multipletests

    pv_array = np.asarray(pv, dtype=float)
    res = np.ones(len(pv_array))

    if not adjust:
        pp = pv_array.copy()
        qpv = np.quantile(pp, umbral)
    else:
        _, pp, _, _ = multipletests(pv_array, method=method)
        qpv = umbral

    mask = pp >= qpv
    res[mask] = (np.log(pp[mask]) / np.log(qpv)) ** alpha

    if isinstance(pv, pd.Series):
        return pd.Series(res, index=pv.index)
    return res


# ---------------------------------------------------------------------------
# nds  (Network Druggability Score)
# ---------------------------------------------------------------------------

def nds(sta, seed, lambda_=0.6, hybrid_gamma=1, cat_rs=None, beta=None):
    """
    Propagate resources from seed targets through annotation categories and
    return a Network Druggability Score for every target.

    Parameters
    ----------
    sta         : DataFrame  species-target-annotation table  [target_id, ann]
    seed        : DataFrame  seed table                       [target_id, w]
    lambda_     : float | None
                  propagation exponent.  Pass None (or np.nan) for hybrid mode,
                  where lambda is computed per target as (k/max_k)^hybrid_gamma.
    hybrid_gamma: float  exponent for dynamic hybridisation (used when lambda_ is None/nan)
    cat_rs      : DataFrame | None  annotation relevance scores  [ann, RS]
    beta        : float | None      exponent for category size correction

    Returns
    -------
    DataFrame  [target_id, score]
    """

    # k_i : number of annotations per target
    kt = sta.groupby("target_id").size().reset_index(name="k")

    use_hybrid = (lambda_ is None) or (isinstance(lambda_, float) and np.isnan(lambda_))
    if use_hybrid:
        k_max = kt["k"].max()
        kt["lambda"] = (kt["k"] / k_max) ** hybrid_gamma
    else:
        kt["lambda"] = lambda_

    # Category sizes and base fac_l = 1 / |category|
    cats = sta.groupby("ann").size().reset_index(name="k")
    cats["fac_l"] = 1.0 / cats["k"]

    # Beta correction on category factor
    if beta is not None:
        cats["fac_l"] = cats["fac_l"] ** beta

    # Relevance score modulation
    if cat_rs is not None:
        cats = cats.merge(cat_rs[["ann", "RS"]], on="ann", how="left")
        cats["fac_l"] = cats["fac_l"] * cats["RS"].fillna(1.0)
        cats = cats.drop(columns=["RS"])

    # fac_seed_j = w_j / k_j ^ lambda_j
    aux = kt.merge(seed, on="target_id", how="inner")
    if use_hybrid:
        aux["fac_seed"] = aux["w"] / (aux["k"] ** aux["lambda"])
    else:
        aux["fac_seed"] = aux["w"] / (aux["k"] ** lambda_)

    # fac_ann_l = sum_j ( fac_seed_j * a_jl )   [sum over seeds annotated to l]
    aux2 = sta.merge(aux[["target_id", "fac_seed"]], on="target_id", how="inner")
    fac_ann = aux2.groupby("ann")["fac_seed"].sum().reset_index(name="fac_ann")

    # f_l = fac_l * fac_ann_l
    f_l = cats.merge(fac_ann, on="ann", how="inner")
    f_l["w"] = f_l["fac_l"] * f_l["fac_ann"]

    # sum_l ( a_il * w_l )  for every target i
    suml_i = sta.merge(f_l[["ann", "w"]], on="ann", how="inner")
    suml_i = suml_i.groupby("target_id")["w"].sum().reset_index(name="suml")

    # score_i = suml_i * k_i ^ (lambda_i - 1)
    score = kt.merge(suml_i, on="target_id", how="inner")
    if use_hybrid:
        score["score"] = score["suml"] * (score["k"] ** (score["lambda"] - 1))
    else:
        score["score"] = score["suml"] * (score["k"] ** (lambda_ - 1))

    score["target_id"] = score["target_id"].astype(str)

    # Add all targets from sta (with score=0 for those not reached by propagation)
    all_targets = pd.DataFrame({"target_id": sta["target_id"].astype(str).unique()})
    score = all_targets.merge(score[["target_id", "score"]], on="target_id", how="left")
    score["score"] = score["score"].fillna(0.0)

    return score.sort_values("score", ascending=False).reset_index(drop=True)
