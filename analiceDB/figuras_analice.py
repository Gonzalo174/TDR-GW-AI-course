"""
figuras_analice.py — las figuras de `analiceDB/`, dibujadas desde sus tablas.

Cada función recibe `t`, el cargador de `comun/figuras.py`, y sólo puede leer
las tablas que declara en `lee`. Ninguna toca la base: las tablas las escribe la
celda de corrida de cada notebook. Se regeneran con

    python comun/figuras.py analiceDB
"""
import numpy as np

import tdr
from figuras import figura

A = "analiceDB"


# --- 01 · dimensiones de la red --------------------------------------------

@figura(A, "01_f01_disponibilidad",
        verifica='Tabla idéntica a la del oráculo v4 (`verificacion/10`); `test_resultados` fija el total de druggables.',
        lee="01_disponibilidad_especie.csv")
def disponibilidad(t, plt):
    """Proteínas vs druggables por especie (log-log)."""
    disp = t("01_disponibilidad_especie.csv")
    d = disp[disp["sp_id"].isin(tdr.ESPECIES_16) & (disp["druggables"] > 0)]
    fig, ax = plt.subplots(figsize=(5.5, 4.5), tight_layout=True)
    ax.scatter(d["proteinas"], d["druggables"], color=tdr.S1, zorder=3)
    for _, r in d.iterrows():
        ax.annotate(tdr.NOMBRE_CORTO.get(r["sp_id"], r["sp_id"]),
                    (r["proteinas"] * 1.04, r["druggables"]), fontsize=7, color=tdr.INK2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("proteínas")
    ax.set_ylabel("druggables")
    ax.set_title("Disponibilidad de evidencia por especie")
    return fig


@figura(A, "01_f02_grado_afiliaciones",
        verifica='`test_resultados.test_la_capa_de_anotaciones_no_se_movio`: categorías y afiliaciones idénticas a las de antes del recorte.',
        lee="01_grado_afiliaciones.csv")
def grado(t, plt):
    """Distribución de grado de las categorías, por fuente (cola pesada: las
    categorías grandes son las que beta > 0 atenúa)."""
    g_all = t("01_grado_afiliaciones.csv")
    fig, ax = plt.subplots(figsize=(5.5, 4), tight_layout=True)
    bins = np.logspace(0, 4, 19)
    for db, etiqueta, color in [("ip", "InterPro", tdr.S1), ("omcl", "OrthoMCL", tdr.S2)]:
        g = g_all.loc[g_all["db"] == db, "grado"]
        ax.hist(g, bins=bins, histtype="step", color=color, label=f"{etiqueta} ({len(g)})")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("proteínas anotadas a la categoría")
    ax.set_ylabel("categorías")
    ax.legend(fontsize=8)
    return fig


@figura(A, "01_f03_tramos_de_grado",
        verifica='La misma tabla que `01_f02`, con el mismo control.',
        lee="01_grado_afiliaciones.csv")
def tramos_de_grado(t, plt):
    """Categorías por tramo de grado y por fuente."""
    grad = t("01_grado_afiliaciones.csv")
    fig, ax = plt.subplots(figsize=(6.5, 3.5), tight_layout=True)
    tab = grad.groupby(["tramo", "db"], observed=True).size().unstack(fill_value=0)
    # el tramo se lee del CSV como texto ("1-5", ">1000"): se ordena por su inicio
    inicio = tab.index.to_series().map(lambda x: int(x.split("-")[0].lstrip(">")))
    tab = tab.loc[inicio.sort_values().index]
    x = np.arange(len(tab))
    ax.bar(x - 0.2, tab.get("ip", 0), width=0.4, color=tdr.S1, label="InterPro")
    ax.bar(x + 0.2, tab.get("omcl", 0), width=0.4, color=tdr.S2, label="OrthoMCL")
    ax.set_xticks(x)
    ax.set_xticklabels(tab.index, fontsize=8)
    ax.set_yscale("log")
    ax.set_xlabel("proteínas anotadas a la categoría")
    ax.set_ylabel("categorías")
    ax.legend(fontsize=8)
    return fig


# --- 02 · capas y conectividad ---------------------------------------------

@figura(A, "02_f01_tamanos_cluster",
        verifica='Sin control externo: describe la capa química tal como quedó tras el recorte, cuyos conteos registra `DB/meta.json`.',
        lee="02_tamanos_cluster.csv")
def tamanos_cluster(t, plt):
    """Compuestos por cluster de fingerprint: sobre qué se corre la capa química."""
    tam = t("02_tamanos_cluster.csv")
    tam = tam[tam["capa"] == "fingerprint"]
    fig, ax = plt.subplots(figsize=(5, 3.5), tight_layout=True)
    ax.hist(tam["tamano"], bins=np.logspace(0, 3, 13), color=tdr.S1)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("compuestos por cluster de fingerprint")
    ax.set_ylabel("clusters")
    ax.set_title("Sobre qué se corre: la capa química")
    return fig


@figura(A, "02_f02_componentes_conexas",
        verifica='Es una verificación del recorte en sí: las dos series tienen que coincidir, porque no debe quedar ninguna componente sin bioactividad positiva.',
        lee="02_componentes_conexas.csv")
def componentes_conexas(t, plt):
    """Tamaño de las componentes conexas, todas vs las que tienen bioactividad
    positiva. Tras el recorte las dos series coinciden."""
    comp = t("02_componentes_conexas.csv")
    fig, ax = plt.subplots(figsize=(5.5, 4), tight_layout=True)
    bins = np.logspace(0, 6, 19)
    ax.hist(comp["total_clusters"], bins=bins, color=tdr.MUTED, label="todas")
    ax.hist(comp.loc[comp["clusters_con_bioact"] > 0, "total_clusters"], bins=bins,
            color=tdr.S1, label="con bioactividad positiva")
    ax.set_xscale("log")
    ax.set_xlabel("clusters en la componente conexa")
    ax.set_ylabel("componentes")
    ax.legend(fontsize=8)
    return fig


@figura(A, "02_f03_conectividad_anotaciones",
        verifica='Sin control: contra el oráculo no hay columnas numéricas comparables (`sin datos` en `10_equivalencia`).',
        lee="02_conectividad_anotaciones.csv")
def conectividad_anotaciones(t, plt):
    """Matriz especie x especie de categorías de afiliación compartidas."""
    matriz = t("02_conectividad_anotaciones.csv", index_col=0)
    fig, ax = plt.subplots(figsize=(6, 5), tight_layout=True)
    im = ax.imshow(np.log10(matriz.values + 1), cmap="Blues")
    etiquetas = [tdr.NOMBRE_CORTO.get(int(s), s) if str(s).isdigit() else s
                 for s in matriz.columns]
    ax.set_xticks(range(len(etiquetas)))
    ax.set_xticklabels(etiquetas, rotation=90, fontsize=7)
    ax.set_yticks(range(len(etiquetas)))
    ax.set_yticklabels(etiquetas, fontsize=7)
    ax.grid(False)
    ax.set_title("Categorías de afiliación compartidas")
    cbar = fig.colorbar(im, ax=ax, label="Cantidad compartida")
    ticks = [0, 1, 10, 100, 1000, 10000]      # la escala es log10(n + 1)
    cbar.set_ticks([np.log10(v + 1) for v in ticks])
    cbar.ax.set_yticklabels([str(v) for v in ticks])
    return fig


# orden alrededor del círculo: las especies de un mismo grupo quedan juntas, y
# el mismo orden en las dos capas deja comparar un grafo con el otro
_ORDEN_GRUPOS = ["Bacterias", "Protozoos", "Amebozoos", "Hongos", "Plantas",
                 "Invertebrados", "Mamiferos"]


def _grafo_conectividad(matriz, color, ax, titulo, que, compartidos):
    """Grafo especie-especie desde una matriz de intersección (análisis de v3).

    La diagonal es el total de cada especie (|A ∩ A| = |A|) y da el área del
    nodo; fuera de la diagonal, lo compartido, que da el ancho y la opacidad de
    la arista. Las dos escalas son lineales y relativas al máximo de la capa.
    Se dibujan todas las aristas con algo en común, en una disposición circular
    fija por grupo taxonómico. (v3 filtraba las aristas por peso relativo
    > 0.4 y usaba kamada-kawai; con la disposición fija el filtro no hace falta,
    y la de v3 cambiaba de un grafo al otro.)
    """
    sp = [int(s) for s in matriz.columns]
    m = matriz.values
    orden = sorted(sp, key=lambda s: (_ORDEN_GRUPOS.index(tdr.SPECIES[s][1]), s))
    ang = np.pi / 2 - 2 * np.pi * np.arange(len(orden)) / len(orden)
    pos = {s: (np.cos(a), np.sin(a)) for s, a in zip(orden, ang)}

    total = {s: m[i, i] for i, s in enumerate(sp)}
    aristas = [(sp[i], sp[j], m[i, j]) for i in range(len(sp))
               for j in range(i + 1, len(sp)) if m[i, j] > 0]
    tot_max = max(total.values())
    com_max = max(c for _, _, c in aristas)
    # de la más fina a la más gruesa: las fuertes quedan arriba
    for u, v, c in sorted(aristas, key=lambda a: a[2]):
        f = c / com_max
        ax.plot(*zip(pos[u], pos[v]), color=color, lw=0.3 + 7 * f, alpha=0.08 + 0.72 * f,
                solid_capstyle="round", zorder=1)
    for s in sp:
        x, y = pos[s]
        ax.scatter(x, y, s=30 + 800 * total[s] / tot_max, color=color,
                   edgecolors="white", linewidths=1, zorder=2)
        # la etiqueta, afuera del círculo y anclada hacia afuera
        ha = "left" if x > 0.3 else "right" if x < -0.3 else "center"
        va = "bottom" if y > 0.5 else "top" if y < -0.5 else "center"
        r = 1.3 if abs(x) < 0.3 else 1.19     # arriba y abajo, escalonadas
        ax.annotate(tdr.NOMBRE_CORTO[s], (r * x, r * y), ha=ha, va=va,
                    fontsize=10.5, color=tdr.INK, style="italic")
    mil = lambda v: f"{int(v):,}".replace(",", "\u2009")
    ax.text(0, -1.55, f"nodo: {que} de la especie (máx. {mil(tot_max)})\n"
            f"arista: {que} {compartidos} (máx. {mil(com_max)})",
            ha="center", va="top", fontsize=9, color=tdr.INK2)
    ax.set_title(titulo)
    ax.set_xlim(-2.0, 2.0)
    ax.set_ylim(-1.9, 1.52)
    ax.set_aspect("equal")
    ax.set_axis_off()


@figura(A, "02_f04_grafo_anotaciones",
        verifica='La misma tabla que `02_f03`, con su mismo alcance: sin control externo (`sin datos` en `10_equivalencia`). La diagonal de la matriz es el total de cada especie.',
        lee="02_conectividad_anotaciones.csv")
def grafo_anotaciones(t, plt):
    """Grafo de especies unidas por categorías de afiliación compartidas: la
    capa de anotaciones, por donde se propaga entre organismos."""
    matriz = t("02_conectividad_anotaciones.csv", index_col=0)
    fig, ax = plt.subplots(figsize=(4.6, 3.9), tight_layout=True)
    _grafo_conectividad(matriz, tdr.S1, ax, "Capa de anotaciones", "categorías",
                        "compartidas")
    return fig


@figura(A, "02_f05_grafo_drogas",
        verifica='La misma tabla que la matriz de drogas de `02`: sin control externo (`sin datos` en `10_equivalencia`). La diagonal de la matriz es el total de cada especie.',
        lee="02_conectividad_drogas.csv")
def grafo_drogas(t, plt):
    """Grafo de especies unidas por compuestos con bioactividad positiva en las
    dos: la capa química."""
    matriz = t("02_conectividad_drogas.csv", index_col=0)
    fig, ax = plt.subplots(figsize=(4.6, 3.9), tight_layout=True)
    _grafo_conectividad(matriz, tdr.S2, ax, "Capa química", "compuestos", "compartidos")
    return fig


# --- 03 · calidad de bioactividades ----------------------------------------

@figura(A, "03_f01_bioactividades_por_tag",
        verifica='Contra el oráculo v4: `cluster_consistent` idéntico; los conteos cambian como explica el recorte.',
        lee="03_bioactividades_por_tag.csv")
def bioactividades_por_tag(t, plt):
    """Registros compuesto-proteína por tag de actividad."""
    tags = t("03_bioactividades_por_tag.csv")
    fig, ax = plt.subplots(figsize=(6, 3.5), tight_layout=True)
    s = tags.groupby("activity_tag")["n"].sum().sort_values()
    ax.barh(range(len(s)), s.values, color=tdr.S1)
    ax.set_yticks(range(len(s)))
    # los tags se muestran con su etiqueta, no con el codigo
    ax.set_yticklabels([tdr.TAG_NOMBRE.get(i, i) for i in s.index], fontsize=8)
    ax.set_xscale("log")
    ax.set_xlabel("registros compuesto-proteína")
    ax.set_title("Composición por tag de actividad")
    return fig


@figura(A, "03_f02_promiscuidad",
        verifica='La lista y la curva vienen de `datos_externos/promiscuidad/derivar.py`, que las reproduce byte a byte con los datos crudos. Contra el oráculo v4: de 25 a 5 promiscuos (compuestos distintos; las tablas tenían filas duplicadas hasta el 2026-09-24), explicado por el recorte. El efecto sobre la base lo mide `03_impacto_en_la_base`.',
        lee="03_curva_filtrado.csv")
def promiscuidad(t, plt):
    """Sensibilidad del criterio de promiscuidad: compuestos y relaciones de
    subestructura que filtraría cada combinación de (MW, N_parentales)."""
    curva = t("03_curva_filtrado.csv")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4), tight_layout=True)
    for umbral, color in zip(sorted(curva["umbral_promiscuidad"].unique()),
                             tdr.ORD + [tdr.S2]):
        c = curva[curva["umbral_promiscuidad"] == umbral]
        ax1.plot(c["mw_max"], c["compuestos_filtrables"], color=color,
                 label=f"N_parentales > {umbral}")
        ax2.plot(c["mw_max"], 100 * c["aristas_filtrables"] / c["aristas_filtrables"].max(),
                 color=color, label=f"N_parentales > {umbral}")
    for ax in (ax1, ax2):
        ax.axvline(tdr.MW_PROMISCUIDAD, color=tdr.MUTED, lw=1, ls="--")
        ax.set_xlabel("umbral de peso molecular (Da)")
        ax.legend(fontsize=8)
    ax1.set_yscale("symlog")
    ax1.set_ylabel("compuestos filtrables")
    ax1.set_title("Compuestos que saca el criterio")
    ax2.set_ylabel("relaciones filtrables (% del máximo)")
    ax2.set_title("Sensibilidad al umbral")
    return fig
