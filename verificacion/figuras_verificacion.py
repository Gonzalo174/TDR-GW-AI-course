"""
figuras_verificacion.py — la figura de la comparación contra `oraculo_v4/`.

Lee `resultados/verificacion_out/10_equivalencia.csv`, que escribe la celda de
corrida de `10_equivalencia.ipynb`. Se regenera con

    python comun/figuras.py verificacion
"""
import numpy as np

import tdr
from figuras import figura

ORDEN = ["idéntico", "explicado por el recorte", "DISCREPANCIA", "sin datos", "sin par nuevo"]
COLOR = {"idéntico": tdr.ST_GOOD, "explicado por el recorte": tdr.S4,
         "DISCREPANCIA": tdr.ST_CRIT, "sin datos": tdr.MUTED, "sin par nuevo": tdr.INK2}


@figura("verificacion", "10_f01_equivalencia",
        verifica='Es la verificación: compara cada tabla contra `oraculo_v4/` con tolerancia relativa 1e-9 (`verificacion/funciones_verificacion.py`).',
        lee="10_equivalencia.csv")
def equivalencia(t, plt):
    """Columnas comparadas contra el oráculo, por análisis y por estado."""
    eq = t("10_equivalencia.csv")
    piv = eq.groupby(["analisis", "estado"]).size().unstack(fill_value=0)
    piv = piv.reindex(columns=[c for c in ORDEN if c in piv.columns])
    fig, ax = plt.subplots(figsize=(7, 3.2), tight_layout=True)
    abajo = np.zeros(len(piv))
    for c in piv.columns:
        ax.barh(range(len(piv)), piv[c], left=abajo, color=COLOR[c], label=c)
        abajo += piv[c].values
    ax.set_yticks(range(len(piv)))
    ax.set_yticklabels(piv.index)
    ax.set_xlabel("columnas comparadas")
    ax.set_title("Equivalencia contra el oráculo, por análisis")
    ax.legend(fontsize=7, loc="lower right")
    return fig
