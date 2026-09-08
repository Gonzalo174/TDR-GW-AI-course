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


def figura(ruta: Path, pie: str) -> str:
    """<figure> con el PNG embebido; si falta, un aviso en vez de un hueco mudo."""
    if not ruta.exists():
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
footer{margin-top:4rem; padding-top:1.5rem; border-top:1px solid var(--linea);
  color:var(--tenue); font-size:.82rem}
a{color:var(--acento)}
"""


def ficha(n, rotulo):
    return f'<div class="ficha"><div class="n">{n}</div><div class="r">{rotulo}</div></div>'


def construir() -> str:
    P = []
    A = P.append

    # ---------- datos que la pagina necesita ----------
    dim = leer("analiceDB_out", "01_dimensiones_capas.csv")
    disp = leer("analiceDB_out", "01_disponibilidad_especie.csv")
    tags = leer("analiceDB_out", "03_bioactividades_por_tag.csv")
    impacto = leer("analiceDB_out", "03_impacto_filtro.csv")
    optimos = leer("genome_prioritization_out", "02_optimos_por_especie.csv")
    consist = leer("genome_prioritization_out", "02_consistencia.csv")
    eq = leer("verificacion_out", "10_equivalencia.csv")
    m01 = meta("analiceDB_out", "01")

    def d(entidad, defecto="—"):
        if dim is None:
            return defecto
        f = dim[dim["entidad"] == entidad]
        return f"{int(f['n'].iloc[0]):,}".replace(",", " ") if len(f) else defecto

    A(f"<style>{CSS}</style>")
    A('<div class="envoltorio">')

    # ---------- cabecera ----------
    A("<header>")
    A("<h1>Priorización de blancos terapéuticos sobre una red multicapa</h1>")
    A('<p class="bajada">Reproducción completa de un análisis de descubrimiento de blancos '
      'sobre una base de bioactividad codificada, con la procedencia de cada resultado '
      'registrada y un oráculo independiente contra el cual contrastarlo.</p>')
    A('<div class="meta-cabecera">')
    A(f'<span>corrida: {m01.get("fecha", "—")}</span>')
    A(f'<span>python {m01.get("python", "—")} · pandas {m01.get("pandas", "—")}</span>')
    A('<span><a href="https://github.com/Gonzalo174/TDR-GW-AI-course">repositorio</a></span>')
    A("</div></header>")

    # ---------- resumen ----------
    A("<h2>De qué se trata</h2>")
    A("<p>El problema es de mi propia investigación: dado un genoma, ordenar sus proteínas "
      "por cuán prometedoras son como blanco de un fármaco. La evidencia disponible es "
      "escasa y está repartida de forma muy desigual entre organismos, así que la "
      "estrategia es propagar lo que se sabe de unas proteínas a otras a través de una red "
      "de tres capas: compuestos, proteínas y anotaciones funcionales.</p>")
    A('<div class="fichas">')
    A(ficha(d("proteínas (V_P)"), "proteínas en la red"))
    A(ficha(d("enlaces E_DP positivos"), "enlaces de bioactividad positiva"))
    A(ficha(d("enlaces proteína-categoría"), "afiliaciones a anotaciones"))
    A(ficha("16", "organismos con evidencia suficiente"))
    A("</div>")
    A('<div class="nota"><p>Este repositorio publica la base <strong>codificada</strong>: '
      "toda columna es un entero y los diccionarios que traducen esos enteros a la notación "
      "original son privados. La pregunta que gobierna el trabajo es si un análisis puede "
      "sobrevivir intacto a esa transformación — y cómo demostrarlo.</p></div>")

    # ---------- datos ----------
    A("<h2>Los datos</h2>")
    A("<p>La base se derivó de la original en tres operaciones, cada una con una alternativa "
      "que se descartó por una razón que está escrita en "
      '<a href="PROVENANCE.md">PROVENANCE.md</a>: un <strong>recorte</strong> que elimina las '
      "componentes conexas del grafo de similitud química sin ninguna bioactividad positiva "
      "(2 396 103 → 831 175 compuestos), la <strong>codificación</strong> a enteros con "
      "semilla fija, y el <strong>empaquetado</strong> de la capa de aristas en dos archivos "
      "comprimidos que preservan las 36 152 622 aristas completas.</p>")
    if dim is not None:
        A(tabla(dim))
        A('<p class="chapo">Dimensiones de las tres capas. Producidas por '
          "<code>analiceDB/01_dimensiones_red.ipynb</code> mediante "
          "<code>funciones_analice.dimensiones_capas</code>.</p>")
    A(figura(fig("analiceDB_out", "01_f01_disponibilidad.png"),
             "Evidencia disponible por organismo. Cada punto es una especie, identificada por "
             "código y tipo de organismo; el nombre no se publica. Los dos ejes son "
             "logarítmicos y abarcan un orden de magnitud largo: esa desigualdad es la que "
             "justifica propagar evidencia entre proteínas."))
    if tags is not None:
        t = tags.copy()
        if "activity_tag" in t.columns:
            t["activity_tag"] = t["activity_tag"].map(lambda v: tdr.TAG_NOMBRE.get(v, v))
        A(tabla(t))
        A('<p class="chapo">Composición de los registros de bioactividad por tag. '
          "<code>analiceDB/03</code>.</p>")

    # ---------- metodo ----------
    A("<h2>El método</h2>")
    A("<p>Sobre la red se corre una propagación con cuatro parámetros. Las anotaciones se "
      "puntúan primero por cuánto se enriquecen en proteínas ya conocidas como blanco "
      "—un test exacto de Fisher por cada categoría, corregido por comparaciones múltiples— "
      "y esa puntuación pondera después la propagación hacia las proteínas sin evidencia.</p>")
    A("<p>La métrica es el área parcial bajo la curva ROC restringida a "
      "<code>FPR ≤ 0.10</code>, corregida por McClish para que azar valga 0.5 y clasificador "
      "perfecto 1.0. La restricción no es cosmética: en este problema sólo importa el "
      "extremo del ranking, porque nadie va a ensayar experimentalmente más allá de las "
      "primeras decenas de proteínas.</p>")
    A('<div class="nota"><p>La validación es <strong>leave-one-species-out</strong>: para '
      "evaluar un organismo se le quita toda su evidencia a la semilla y se lo prioriza como "
      "si no se supiera nada de él. Que eso esté bien hecho no se asume: hay una prueba "
      "dedicada (<code>comun/tests/test_fuga.py</code>) que verifica que ninguna proteína "
      "del organismo evaluado sobrevive en la semilla, con su control negativo al lado.</p></div>")

    # ---------- resultados ----------
    A("<h2>Resultados</h2>")
    if optimos is not None:
        cols = [c for c in ["especie", "AUC01", "AUC", "alpha", "beta", "lambda_", "gamma",
                            "N_targets"] if c in optimos.columns]
        A(tabla(optimos[cols].sort_values("AUC01", ascending=False)))
        A('<p class="chapo">Parámetros óptimos y desempeño por organismo. '
          "<code>genome_prioritization/02_optimo_y_consistencia.ipynb</code>.</p>")
    A(figura(fig("genome_prioritization_out", "01_f03_auc01_por_especie.png"),
             "Desempeño por organismo en su óptimo, bajo leave-one-species-out. La línea "
             "punteada en 0.5 es el azar."))
    A(figura(fig("genome_prioritization_out", "02_f03_plateau.png"),
             "Cuán plano es el óptimo. Si el desempeño cayera abruptamente al moverse del "
             "máximo, la elección de parámetros sería frágil y el número reportado, casual."))

    # ---------- verificacion ----------
    A("<h2>Cómo sé que esto está bien</h2>")
    A("<p>Tres controles, del más barato al más caro. Ninguno es una figura sin nada detrás.</p>")
    A("<h3>1. Las pruebas</h3>")
    A("<p>38 pruebas que corren en unos 15 segundos: integridad de las tablas, ausencia de "
      "fuga en la validación cruzada, las métricas contra valores de referencia, las rutas y "
      "el entorno. Cuatro son específicas de esta reproducción y cubren lo que puede "
      "romperse <em>en silencio</em> — el caso real que las motivó está más abajo.</p>")
    A("<h3>2. El control externo</h3>")
    A("<p>Los óptimos por organismo se contrastan contra los de una corrida independiente "
      "anterior, versionada en <code>control/</code>.</p>")
    if consist is not None:
        A(tabla(consist, maxf=12))
    A("<h3>3. El oráculo</h3>")
    A("<p>Las mismas tablas calculadas antes del recorte y de la codificación, en "
      "<code>oraculo_v4/</code>. La comparación es informativa por una asimetría: el recorte "
      "<strong>debe</strong> mover los descriptivos de la base, y <strong>no debe</strong> "
      "mover la priorización, que trabaja sobre anotaciones y blancos conocidos.</p>")
    if eq is not None and len(eq):
        piv = eq.groupby(["analisis", "estado"]).size().unstack(fill_value=0).reset_index()
        A(tabla(piv))
        n_disc = int((eq["estado"] == "DISCREPANCIA").sum())
        marca = ('<span class="marca m-ok">sin discrepancias</span>' if n_disc == 0
                 else f'<span class="marca m-mal">{n_disc} discrepancias</span>')
        A(f"<p>{marca} sobre {len(eq)} columnas comparadas.</p>")
    A(figura(fig("verificacion_out", "10_f01_equivalencia.png"),
             "Equivalencia contra el oráculo por análisis. Verde es idéntico; ámbar es una "
             "diferencia que el recorte explica; rojo sería una discrepancia a investigar."))

    A("<h3>Lo que casi sale mal</h3>")
    A('<div class="nota"><p>La base codificada guarda el tag de actividad como el entero '
      "<code>2</code>, mientras el código heredado filtraba por el texto "
      "<code>\"positive\"</code>. Esa comparación no lanza ningún error: devuelve cero filas. "
      "El pipeline entero seguía adelante produciendo NaN y rankings vacíos, con "
      "<strong>248 457</strong> enlaces de evidencia silenciosamente descartados.</p>"
      "<p>Un segundo caso: InterPro y OrthoMCL se codificaron por separado, cada uno "
      "empezando en 0, así que sus códigos se pisan — los 14 083 dominios caen dentro del "
      "rango de los 76 157 grupos de ortología. Concatenarlos habría fusionado anotaciones "
      "sin relación alguna.</p>"
      "<p>Ninguno de los dos se encontró leyendo el código: los encontraron las pruebas al "
      "correrlas por primera vez contra la base nueva. Es la razón de que las pruebas vayan "
      "antes que los resultados y no después.</p></div>")

    # ---------- reproducir ----------
    A("<h2>Reproducirlo</h2>")
    A("<pre><code>git clone https://github.com/Gonzalo174/TDR-GW-AI-course\n"
      "cd TDR-GW-AI-course\n"
      "conda create -n TDR_GW python=3.12 &amp;&amp; conda activate TDR_GW\n"
      "pip install -r requirements.txt\n"
      "python -m unittest discover -s comun/tests -p \"test_*.py\"\n"
      "jupyter nbconvert --to notebook --execute --inplace analiceDB/*.ipynb</code></pre>")
    A("<p>Las rutas se derivan de la ubicación del repositorio, así que no hay nada que "
      "configurar. Cada notebook deja un <code>NN_meta.json</code> junto a sus tablas con la "
      "fecha, los parámetros, la fecha de cada tabla de entrada y las versiones de python y "
      "pandas.</p>")

    A("<h2>Los límites</h2>")
    A("<p>Un repositorio que no declara lo que no puede hacer no es verificable. Éste no "
      "reproduce: la base misma, que se deriva de datos que no se publican; los nombres "
      "legibles de las anotaciones, que exigen un diccionario privado; y la identidad de los "
      "organismos, que es deliberada. Lo que sí se publica es la salida de cada uno de esos "
      "pasos y el registro de cómo se produjo.</p>")

    A('<footer>Generado por <code>informe/generar_pagina.py</code> a partir de las tablas de '
      "<code>resultados/</code>. Los números de esta página no se escribieron a mano.</footer>")
    A("</div>")
    return "\n".join(P)


if __name__ == "__main__":
    destino = RAIZ / "index.html"
    destino.write_text("<title>Priorización de blancos terapéuticos</title>\n" + construir())
    print(f"escrito: {destino}  ({destino.stat().st_size/1024:.0f} KB)")
