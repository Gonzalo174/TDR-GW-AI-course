#!/usr/bin/env python3
"""Genera `index.html` a partir de las salidas de los notebooks.

La pagina no se escribe a mano: se arma leyendo las tablas de `resultados/` y
embebiendo las figuras que produjeron los notebooks. Si un numero cambia porque
se volvio a correr un analisis, la pagina cambia con el; no hay forma de que el
texto diga una cosa y la tabla otra.

    python3 informe/generar_pagina.py        # -> index.html

Las figuras van embebidas en base64 para que la pagina sea un unico archivo que
se pueda linkear desde cualquier lado.
"""
from __future__ import annotations

import base64
import html
import json
import sys
from pathlib import Path

import pandas as pd

RAIZ = Path(__file__).resolve().parents[1]
RES = RAIZ / "resultados"
sys.path.insert(0, str(RAIZ / "comun"))
import tdr  # noqa: E402


# --------------------------------------------------------------------------
# Piezas
# --------------------------------------------------------------------------

def tabla(df: pd.DataFrame, clases="", maxf=None) -> str:
    """DataFrame -> <table>. Sin indice, con los numeros ya formateados."""
    if maxf:
        df = df.head(maxf)
    df = df.copy()
    for c in df.columns:
        if pd.api.types.is_float_dtype(df[c]):
            df[c] = df[c].map(lambda v: "" if pd.isna(v) else f"{v:,.4g}")
        elif pd.api.types.is_integer_dtype(df[c]):
            df[c] = df[c].map(lambda v: f"{v:,}".replace(",", " "))
    th = "".join(f"<th>{html.escape(str(c))}</th>" for c in df.columns)
    tr = "".join("<tr>" + "".join(f"<td>{html.escape(str(v))}</td>" for v in fila) + "</tr>"
                 for fila in df.itertuples(index=False))
    return f'<div class="tabla-scroll"><table class="{clases}"><thead><tr>{th}</tr></thead><tbody>{tr}</tbody></table></div>'


FALTANTES = []   # figuras que no estaban: si hay alguna, la pagina no se escribe


def figura(ruta: Path, pie: str) -> str:
    """<figure> con el PNG embebido; si falta, un aviso en vez de un hueco mudo."""
    if not ruta.exists():
        FALTANTES.append(ruta.name)
        return f'<figure class="falta"><div class="aviso">falta {html.escape(ruta.name)}</div><figcaption>{pie}</figcaption></figure>'
    b64 = base64.b64encode(ruta.read_bytes()).decode()
    return (f'<figure><img src="data:image/png;base64,{b64}" alt="{html.escape(pie[:80])}">'
            f'<figcaption>{pie}</figcaption></figure>')


def leer(sub: str, nombre: str):
    p = RES / sub / nombre
    return pd.read_csv(p) if p.exists() else None


def meta(sub: str, nb: str):
    p = RES / sub / f"{nb}_meta.json"
    return json.loads(p.read_text()) if p.exists() else {}


def fig(sub: str, nombre: str) -> Path:
    return RES / sub / "figuras" / nombre


CSS = """
:root{
  --tinta:#12100c; --tinta2:#4a463d; --tenue:#8a857a;
  --fondo:#fbfaf7; --panel:#ffffff; --linea:#e4e1d8;
  --acento:#2a78d6; --ok:#0ca30c; --alerta:#d03b3b; --ambar:#eda100;
}
@media (prefers-color-scheme: dark){
  :root:not([data-theme="light"]){
    --tinta:#eceae4; --tinta2:#b6b1a6; --tenue:#8a857a;
    --fondo:#131211; --panel:#1b1a18; --linea:#302e2a;
  }
}
:root[data-theme="dark"]{
  --tinta:#eceae4; --tinta2:#b6b1a6; --tenue:#8a857a;
  --fondo:#131211; --panel:#1b1a18; --linea:#302e2a;
}
*{box-sizing:border-box}
body{background:var(--fondo); color:var(--tinta);
  font:16px/1.65 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif;
  margin:0; padding:0 1.2rem 5rem;}
.envoltorio{max-width:60rem; margin:0 auto}
header{padding:3.5rem 0 2rem; border-bottom:1px solid var(--linea); margin-bottom:2.5rem}
h1{font-size:2rem; line-height:1.2; margin:0 0 .6rem; letter-spacing:-.02em}
.bajada{color:var(--tinta2); font-size:1.05rem; margin:0 0 1.2rem; max-width:44rem}
.meta-cabecera{color:var(--tenue); font-size:.82rem; display:flex; flex-wrap:wrap; gap:.4rem 1.2rem}
h2{font-size:1.3rem; margin:3rem 0 .4rem; letter-spacing:-.01em}
h2::after{content:""; display:block; width:2.5rem; height:2px; background:var(--acento); margin-top:.5rem}
h3{font-size:1rem; margin:2rem 0 .3rem; color:var(--tinta)}
p{margin:.7rem 0; max-width:44rem}
.chapo{color:var(--tinta2)}
code{font-family:ui-monospace,SFMono-Regular,Menlo,monospace; font-size:.86em;
  background:var(--panel); border:1px solid var(--linea); border-radius:4px; padding:.05em .35em}
pre{background:var(--panel); border:1px solid var(--linea); border-radius:8px;
  padding:.9rem 1rem; overflow-x:auto; font-size:.84rem; line-height:1.5}
pre code{background:none; border:none; padding:0}
.tabla-scroll{overflow-x:auto; margin:1rem 0; border:1px solid var(--linea); border-radius:8px; background:var(--panel)}
table{border-collapse:collapse; width:100%; font-size:.85rem}
th,td{text-align:left; padding:.45rem .7rem; border-bottom:1px solid var(--linea); white-space:nowrap}
th{background:color-mix(in srgb, var(--linea) 35%, transparent); font-weight:600; color:var(--tinta2)}
tbody tr:last-child td{border-bottom:none}
td:not(:first-child){font-variant-numeric:tabular-nums}
figure{margin:1.6rem 0; background:var(--panel); border:1px solid var(--linea);
  border-radius:10px; padding:1rem; overflow:hidden}
figure img{display:block; width:100%; height:auto; border-radius:4px}
figcaption{color:var(--tinta2); font-size:.84rem; margin-top:.7rem; line-height:1.5}
.aviso{color:var(--alerta); font-size:.85rem; padding:1.5rem; text-align:center;
  border:1px dashed var(--alerta); border-radius:6px}
.fichas{display:grid; grid-template-columns:repeat(auto-fit,minmax(11rem,1fr)); gap:.8rem; margin:1.4rem 0}
.ficha{background:var(--panel); border:1px solid var(--linea); border-radius:10px; padding:.9rem 1rem}
.ficha .n{font-size:1.5rem; font-weight:650; letter-spacing:-.02em; font-variant-numeric:tabular-nums}
.ficha .r{color:var(--tenue); font-size:.78rem; margin-top:.15rem; line-height:1.35}
.marca{display:inline-block; font-size:.72rem; font-weight:600; padding:.12rem .5rem;
  border-radius:99px; border:1px solid currentColor}
.m-ok{color:var(--ok)} .m-cambia{color:var(--ambar)} .m-mal{color:var(--alerta)}
.nota{border-left:3px solid var(--acento); background:var(--panel);
  padding:.7rem 1rem; margin:1.2rem 0; border-radius:0 8px 8px 0}
.nota p{margin:.3rem 0}
.par{display:grid; grid-template-columns:repeat(auto-fit,minmax(20rem,1fr)); gap:0 1rem}
.par figure{margin:1rem 0}
.diapo{background:var(--panel); border:1px solid var(--linea); border-radius:12px;
  padding:1.2rem 1.4rem; display:grid; grid-template-columns:repeat(auto-fit,minmax(22rem,1fr));
  gap:0 1.4rem; align-items:center; margin:1.2rem 0}
.diapo figure{border:none; padding:0; margin:.6rem 0}
.diapo ul{margin:.4rem 0; padding-left:1.1rem; font-size:.93rem}
.diapo li{margin:.35rem 0}
footer{margin-top:4rem; padding-top:1.5rem; border-top:1px solid var(--linea);
  color:var(--tenue); font-size:.82rem}
a{color:var(--acento)}
"""


def ficha(n, rotulo):
    return f'<div class="ficha"><div class="n">{n}</div><div class="r">{rotulo}</div></div>'


def numeros() -> dict:
    """Las mismas cifras que el informe (`informe/numeros.py`), sin el marcado
    de LaTeX: una sola fuente para los dos entregables."""
    sys.path.insert(0, str(RAIZ / "informe"))
    import numeros as nu
    return {k: str(v).replace("\\,", "\u202f") for k, v in nu.construir().items()}


def construir() -> str:
    P = []
    A = P.append
    n = numeros()

    optimos = leer("genome_prioritization_out", "02_optimos_por_especie.csv")
    m01 = meta("analiceDB_out", "01")

    A(f"<style>{CSS}</style>")
    A('<div class="envoltorio">')

    # ---------- cabecera ----------
    A("<header>")
    A("<h1>Priorización de blancos terapéuticos sobre una red multicapa</h1>")
    A('<p class="bajada">Propagar la poca evidencia experimental que hay entre organismos '
      "para ordenar un genoma completo por cuán prometedoras son sus proteínas como blanco "
      "de un fármaco, y usar la misma red para proponer el blanco de compuestos que no "
      "tienen ninguno.</p>")
    A('<div class="meta-cabecera">')
    A(f'<span>corrida: {m01.get("fecha", "—")}</span>')
    A(f'<span>python {m01.get("python", "—")} · pandas {m01.get("pandas", "—")}</span>')
    A('<span><a href="https://github.com/Gonzalo174/TDR-GW-AI-course">repositorio</a></span>')
    A('<span><a href="informe/informe.pdf">informe (PDF)</a></span>')
    A('<span><a href="informe/presentacion.pdf">presentación (PDF)</a></span>')
    A("</div></header>")

    # ---------- resumen ----------
    A('<p class="chapo" lang="en"><strong>In English.</strong> Given a genome, rank its '
      "proteins as drug-target candidates by propagating scarce experimental evidence across a "
      "three-layer network (compounds, proteins, functional annotations). Evaluated leaving out "
      f'all evidence of the organism being ranked, it beats chance in {n["NSobreAzar"]} of '
      f'{n["NEspeciesOpt"]} organisms (median AUC01 {n["AUCMediana"]}). Run on an '
      "integer-coded database; every number on this page is read from the result tables, and "
      "every figure is redrawn from them in about 40 s.</p>")
    A("<h2>En una pantalla</h2>")
    A('<div class="fichas">')
    A(ficha(f'{n["NSobreAzar"]} de {n["NEspeciesOpt"]}', "organismos priorizados mejor que el "
            "azar, sin usar nada de su propia evidencia"))
    A(ficha(n["AUCMediana"], f'AUC01 mediana (rango {n["AUCMin"]}–{n["AUCMax"]})'))
    A(ficha(f'{n["PctRecInformativa"]} %', "de los blancos recuperados cuando la semilla "
            "comparte una anotación con el blanco"))
    A(ficha(f'{n["PctSinSemilla"]} %', "de los compuestos sin ningún vecino químico con "
            "blanco: el cuello de botella"))
    A("</div>")
    A("<p>El problema es de mi propia investigación: dado un genoma, ordenar sus proteínas "
      "por cuán prometedoras son como blanco de un fármaco. La evidencia es escasa y está "
      "repartida de forma muy desigual entre organismos, así que la estrategia es propagar "
      "lo que se sabe sobre una red de tres capas —compuestos, proteínas y anotaciones "
      "funcionales— desde las proteínas con evidencia hacia las que no la tienen.</p>")

    # ---------- datos y metodo ----------
    A("<h2>Datos y método</h2>")
    A(f'<p>La base tiene {n["NProteinas"]} proteínas de {n["NEspecies"]} organismos con '
      f'evidencia suficiente, {n["NPositivos"]} enlaces de bioactividad positiva y '
      f'{n["NAfiliaciones"]} afiliaciones a {n["NInterPro"]} dominios InterPro y '
      f'{n["NOrthoMCL"]} grupos de ortología. La evidencia por organismo va de '
      f'{n["DrugMin"]} a {n["DrugMax"]} proteínas con actividad conocida: esa desigualdad '
      "es la premisa del método. La base publicada está <strong>codificada</strong> (toda "
      "columna es un entero; los diccionarios son privados) y <strong>recortada</strong> a "
      "las componentes químicas con alguna bioactividad positiva.</p>")
    A(figura(fig("analiceDB_out", "01_f01_disponibilidad.png"),
             "Evidencia disponible por organismo, en escala logarítmica en los dos ejes. Cada "
             "punto es una especie, identificada por código y tipo de organismo."))
    A("<p>Las anotaciones se puntúan por cuánto se enriquecen en blancos conocidos (Fisher "
      "unilateral, corregido por Benjamini–Hochberg) y esa puntuación pondera la propagación, "
      "que tiene cuatro parámetros: α (peso de la puntuación), β (penalización de las "
      "categorías grandes), λ y γ (cómo se combinan las categorías de una proteína). La "
      "métrica es la <strong>AUC01</strong>: el área bajo la ROC hasta <code>FPR ≤ 0.10</code>, "
      "corregida por McClish para que el azar valga 0.5. Sólo importa el extremo del "
      "ranking.</p>")
    A('<div class="nota"><p>La validación es <strong>leave-one-species-out</strong>: para '
      "evaluar un organismo se le quita toda su evidencia a la semilla y se lo prioriza como "
      "si no se supiera nada de él. Una prueba dedicada "
      "(<code>comun/tests/test_fuga.py</code>) verifica que ninguna proteína del organismo "
      "evaluado sobrevive en la semilla.</p></div>")

    # ---------- priorizacion ----------
    A("<h2>1 · Priorización de genoma completo</h2>")
    A(f'<p>Para cada organismo se barrieron las 120 combinaciones válidas de parámetros. En su '
      f'óptimo la AUC01 mediana es <strong>{n["AUCMediana"]}</strong> y los '
      f'{n["NSobreAzar"]} organismos superan el azar. El orden no es arbitrario: los '
      f'{n["NProtozoosArriba"]} primeros puestos son los {n["NProtozoos"]} protozoos '
      "parásitos.</p>")
    A(figura(fig("genome_prioritization_out", "01_f03_auc01_por_especie.png"),
             "AUC01 en el óptimo de cada organismo, bajo leave-one-species-out. La línea "
             "punteada es el azar."))
    if optimos is not None:
        o = optimos.merge(tdr.META[["sp", "grupo", "parasito"]], left_on="especie",
                          right_on="sp", how="left")
        o["especie"] = o["especie"].map(lambda s: f"sp{int(s):02d}")
        o["lambda_"] = o["lambda_"].map(lambda v: "híbrido" if pd.isna(v) else f"{v:g}")
        o["parasito"] = o["parasito"].map({True: "sí", False: "no"})
        cols = ["especie", "grupo", "parasito", "AUC01", "AUC", "alpha", "beta", "lambda_", "gamma"]
        A(tabla(o[cols].sort_values("AUC01", ascending=False)))
        A('<p class="chapo">Óptimo por organismo. '
          "<code>genome_prioritization/02_optimo_y_consistencia.ipynb</code>.</p>")
    A(figura(fig("genome_prioritization_out", "02_f01_roc_26.png"),
             f'El organismo foco, un protozoo parásito (sp26), en su óptimo: AUC global '
             f'{n["AUCGlobalFoco"]}, AUC01 {n["AUCFoco"]}. A la derecha, el puntaje de '
             "propagación de las proteínas druggables frente al resto."))
    A("<h3>Qué decide el resultado</h3>")
    A(f'<p>De los cuatro parámetros sólo <strong>β</strong> mueve la AUC01 de forma '
      f'consistente: penalizar las categorías grandes mejora la priorización, y el óptimo tiene '
      f'β &gt; 0 en {n["NBetaPositivo"]} de {n["NEspeciesOpt"]} organismos. Entre las diez '
      f'mejores combinaciones de cada organismo, β = 1 aparece {n["EnrBetaUno"]} veces más de '
      f'lo que daría el azar y β = 0 apenas {n["EnrBetaCero"]}.</p>')
    A(figura(fig("genome_prioritization_out", "02_f02_grilla_beta.png"),
             "AUC01 de las 120 combinaciones por organismo, agrupadas por el valor de cada "
             "parámetro. En el panel de β, 0 es la red G'r y &gt;0 la G'rk."))
    A("<h3>Qué tan firme es el óptimo</h3>")
    A(f'<p>En todos los organismos las veinte mejores combinaciones quedan a menos del '
      f'{n["CaidaTopVeinteMax"]} % del máximo: el óptimo es una meseta, no un pico casual. Y '
      f'el perfil de la grilla se parece entre organismos (Spearman mediana '
      f'{n["SpearmanMediana"]}, positiva en el {n["PctSpearmanPositivo"]} % de los pares): '
      "los parámetros que funcionan en uno funcionan en los demás.</p>")
    A('<div class="par">')
    A(figura(fig("genome_prioritization_out", "02_f03_plateau.png"),
             "Caída de la AUC01 entre el óptimo y la K-ésima mejor combinación."))
    A(figura(fig("genome_prioritization_out", "02_f04_consistencia_especies.png"),
             "Concordancia entre organismos sobre el perfil de la grilla."))
    A("</div>")

    # ---------- huerfanas ----------
    A("<h2>2 · Desorfanización de compuestos</h2>")
    A(f'<p>La segunda pregunta invierte la dirección: dado un compuesto sin blanco conocido, '
      f'proponer uno. Para medir el acierto se usan <strong>pseudohuérfanas</strong>, compuestos '
      f'con exactamente un blanco conocido a los que se les borra toda su bioactividad; la '
      f'semilla se arma desde su vecindario químico y se mide la posición relativa del blanco '
      f'verdadero (<em>frank</em>, recuperado si &lt; 0.1). Se evaluaron {n["NPseudo"]}.</p>')
    A(f'<p>El blanco se recupera en el <strong>{n["PctRecuperadas"]} %</strong> de los casos, y la '
      f'causa está medida: la recuperación la decide la <em>clase de semilla</em>, no el modelo. '
      f'El {n["PctSinSemilla"]} % no tiene ningún vecino químico con blanco (semilla nula). '
      f'Con semilla pero sin ninguna anotación en común con el blanco '
      f'({n["NNoInformativa"]} casos), la recuperación es del {n["PctRecNoInformativa"]} %; con '
      f'al menos una ({n["NInformativa"]}), del <strong>{n["PctRecInformativa"]} %</strong>.</p>')
    A('<div class="par">')
    A(figura(fig("huerfanas_out", "01_f04_recuperacion_por_semilla.png"),
             "Recuperación por clase de semilla."))
    A(figura(fig("huerfanas_out", "01_f03_semilla_por_especie.png"),
             "Composición de la semilla por organismo: la fracción informativa (verde) es el "
             "techo del método."))
    A("</div>")
    A(f'<p>La fracción informativa es un techo que ninguna reponderación de la red puede '
      f'superar; por organismo va del {n["TechoMin"]} % al {n["TechoMax"]} %. Tres palancas para '
      f'ampliar la semilla (subir el umbral de similitud, sumar KEGG, usar la actividad '
      f'fenotípica) rescataron {n["NRescatadas"]} de {n["NMuestraPalancas"]}; las dos primeras '
      "no podían hacerlo en esta base.</p>")
    A("<h3>Cuando hay semilla informativa</h3>")
    A(f'<p>El {n["PctRecInformativa"]} % dice poco: con unas doce mil proteínas por organismo, '
      f'frank &lt; 0.1 es quedar entre las primeras ~{n["FrankCorteRank"]}. Mirando la posición '
      f'absoluta, el blanco queda <strong>primero en el {n["PctInfTopUno"]} %</strong> de los '
      f'casos, en el top-10 en el {n["PctInfTopDiez"]} % y en el top-100 en el '
      f'{n["PctInfTopCien"]} % (el azar pondría en el top-10 al {n["AzarTopDiez"]} %). Son dos '
      f'poblaciones: en la inferencia <em>directa</em> ({n["NInfDirecta"]}; el blanco ya está en '
      f'la semilla) la posición mediana es {n["RSSMedInfDirecta"]} y el top-10 llega al '
      f'{n["PctTopDiezDirecta"]} %; en la <em>indirecta</em> ({n["NInfIndirecta"]}; sólo por '
      f'anotaciones) la mediana es {n["RSSMedInfIndirecta"]} y el top-10, '
      f'{n["PctTopDiezIndirecta"]} %. Las semillas grandes diluyen la señal (ρ = '
      f'{n["RhoSemillaRSS"]}, positiva en {n["NOrgRhoPos"]} de {n["NOrgRho"]} organismos con '
      f'casos suficientes). En el {n["PctEmpates"]} % el blanco empata en puntaje con otras '
      "proteínas.</p>")
    A(figura(fig("huerfanas_out", "03_f05_semilla_informativa.png"),
             "(a) Fracción con el blanco en el top-k de su organismo, frente al azar; la "
             "punteada es frank = 0.1. (b) frank en escala log. (c) Tamaño de la semilla contra "
             "posición del blanco. (d) Posición por organismo; la barra es la mediana."))
    A("<h3>Hasta dónde confiar en el ranking</h3>")
    A(f'<p>Sobre el ranking global, la recuperación acumulada tiene pendiente alta en las '
      f'primeras posiciones y cae al ruido de fondo en <strong>r*G = {n["RGEstrella"]}</strong>: '
      f'una propuesta más abajo no se distingue del azar. La inferencia <em>directa</em> (el '
      f'blanco lo trae un vecino químico) es el {n["PctDirecta"]} % de los casos y pone el blanco '
      f'en la posición global mediana {n["RGDirecta"]}; la <em>indirecta</em> (sólo por '
      f'anotaciones), el {n["PctIndirecta"]} %, lo pone en la {n["RGIndirecta"]}.</p>')
    A(figura(fig("huerfanas_out", "03_f02_recuperacion.png"),
             "Arriba: compuestos cuyo blanco cae antes de la posición l del ranking global. "
             "Abajo: su derivada suavizada y el umbral λ∞ + 3σ que define r*G."))
    A(f'<p><strong>Aplicación al organismo foco.</strong> De {n["EmbudoActivos"]} compuestos con '
      f'actividad fenotípica contra sp26, {n["EmbudoHuerfanos"]} no tienen blanco proteico y '
      f'{n["EmbudoTratables"]} tienen vecinos químicos con blanco; ninguno recibe una propuesta '
      f'por encima de r*G ({n["NSugerencias"]} sugerencias). Con la base recortada el embudo es '
      "demasiado angosto, y el método lo dice en vez de forzar una respuesta.</p>")

    # ---------- contraste con v4: un solo panel ----------
    A("<h2>3 · Contraste con la versión anterior (v4)</h2>")
    A('<div class="diapo">')
    A('<div class="diapo-texto">')
    A("<p>Este repositorio porta el análisis v4 a la base codificada y recortada. Las mismas "
      "tablas calculadas antes de los dos cambios quedaron como oráculo, con un único criterio: "
      "<strong>si el análisis toca la capa química, el recorte lo mueve; si no, no</strong>. "
      f'Sobre {n["NComparadas"]} columnas: ')
    marca = ('<span class="marca m-ok">sin discrepancias</span>' if n["NDiscrepancias"] == "0"
             else f'<span class="marca m-mal">{n["NDiscrepancias"]} discrepancias</span>')
    A(f"{marca}</p><ul>")
    A(f'<li><strong>Priorización</strong> (sin química): {n["NIdenticasGP"]} de '
      f'{n["NComparadasGP"]} columnas idénticas; el control independiente coincide en '
      f'{n["ControlCoinciden"]} de {n["NEspeciesOpt"]} organismos, Δ AUC01 máx. = '
      f'{n["ControlDeltaMax"]}.</li>')
    A(f'<li><strong>Desorfanización</strong>: {n["NCambianHU"]} de {n["NComparadasHU"]} columnas '
      f'cambian como predice el recorte. Menos compuestos promiscuos filtrados '
      f'({n["PromiscuosVcuatro"]} → {n["NPromiscuos"]}), más vecinos: semilla nula '
      f'{n["PctSinSemillaVcuatro"]} → {n["PctSinSemilla"]} %, recuperación '
      f'{n["PctRecuperadasVcuatro"]} → {n["PctRecuperadas"]} %, r*G '
      f'{n["RGEstrellaVcuatro"]} → {n["RGEstrella"]}.</li>')
    A(f'<li><strong>Descriptivos</strong>: {n["NCambianAN"]} de {n["NComparadasAN"]} cambian, las '
      "que miden componentes químicas.</li></ul>")
    A('<p class="chapo">El port destapó cuatro filtros heredados que comparaban un entero '
      "codificado contra el texto que solía contener; ninguno lanzaba una excepción, y los "
      "encontraron las pruebas y esta comparación. Detalle en "
      '<a href="PROVENANCE.md">PROVENANCE.md</a>.</p>')
    A("</div>")
    A(figura(fig("verificacion_out", "10_f01_equivalencia.png"),
             "Columnas comparadas contra el oráculo, por análisis. Verde: idéntica; ámbar: "
             "cambio explicado por el recorte."))
    A("</div>")

    # ---------- reproducir ----------
    A("<h2>Verificarlo y reproducirlo</h2>")
    A(f'<p>{n["NPruebas"]} pruebas corren en menos de medio minuto: integridad de las tablas, '
      "ausencia de fuga, métricas contra valores de referencia, rutas, entorno, y que toda "
      "figura se dibuje sin tocar la base. Las tablas de la corrida están versionadas, así que "
      f'las {n["NFiguras"]} figuras y todas las cifras de esta página <strong>se reproducen '
      "desde las tablas pregeneradas en alrededor de un minuto</strong>, sin volver a correr "
      "ningún modelo. Rehacer las tablas desde la base lleva unos 45 minutos con 20 "
      "procesos.</p>")
    A("<pre><code>conda create -n TDR_GW python=3.12 &amp;&amp; conda activate TDR_GW\n"
      "pip install -r requirements.txt\n"
      "make verificar     # pruebas + figuras desde las tablas + cifras del informe   ~1 min\n"
      "make corrida       # los 10 notebooks, desde la base                          ~45 min\n"
      "make entregables   # figuras, informe, presentación y esta página         ~1,5 min</code></pre>")
    A("<p>Cada notebook deja un <code>NN_meta.json</code> junto a sus tablas con la fecha, los "
      "parámetros, la fecha de cada tabla de entrada y las versiones de python y pandas. Qué "
      "tablas lee cada figura y qué la verifica está en <code>FIGURAS.md</code>. El repositorio "
      "no reproduce la base misma ni el filtro de promiscuidad de <code>analiceDB/03</code> "
      "(los dos necesitan datos privados; su salida está versionada), ni los nombres legibles "
      "de las anotaciones ni la identidad de los organismos, que es deliberada.</p>")

    A('<footer>Generado por <code>informe/generar_pagina.py</code> a partir de las tablas de '
      "<code>resultados/</code>. Los números de esta página no se escribieron a mano.</footer>")
    A("</div>")
    return "\n".join(P)


if __name__ == "__main__":
    pagina = construir()
    # Como numeros.py: si falta algo no se pisa la pagina versionada.
    sin_dato = pagina.count("??")
    if (FALTANTES or sin_dato) and "--forzar" not in sys.argv:
        print(f"NO se escribe index.html: faltan {len(FALTANTES)} figuras "
              f"({', '.join(FALTANTES[:5])}) y {sin_dato} cifras.\n"
              "Las figuras se dibujan desde las tablas con `python comun/figuras.py`. "
              "Para escribirla igual: --forzar")
        sys.exit(1)
    destino = RAIZ / "index.html"
    destino.write_text("<title>Priorización de blancos terapéuticos</title>\n" + pagina)
    print(f"escrito: {destino}  ({destino.stat().st_size/1024:.0f} KB)")
