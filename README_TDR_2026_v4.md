# TDR_2026_v4 — organización propuesta

Documento de trabajo escrito en `TDR_2026_v3/` para usarse como `README.md` de
`TDR_2026_v4/`. Recoge la estructura pedida (`createDB`, `analiceDB`,
`genome_prioritization`, `huerfanas`), el formato de notebook y el desdoblamiento
de funciones en principales/auxiliares, más los cambios que conviene hacerle al
enfoque antes de empezar a mover archivos.

**Convención de datos.** Las tablas de entrada se leen de `gon3/DB/` (la base no
se regenera en v4). **Todas las salidas escriben en el symlink `gon4`** —nunca en
`gon3/` ni en el home— dentro de una carpeta que repite el nombre del análisis
con el sufijo `_out`: `genome_prioritization/` escribe en
`gon4/genome_prioritization_out/`.

**Convención de numeración.** Las carpetas **no** llevan número (así son
importables desde Python); los notebooks dentro de cada una **sí**, y **cada
archivo que un notebook genera arranca con su mismo número**: el `01` produce
`01_*.csv`, `01_meta.json` y `01_f0N_*.pdf`. Mirando un archivo de salida se sabe
qué notebook hay que abrir para regenerarlo.

Estado al escribirlo (2026-09-04): el barrido de `full_genome_v5.ipynb` terminó
y dejó las 16 especies en `gon3/resultados/genoma_completo_v5/`; esa corrida se
rehace en v4 escribiendo a `gon4/genome_prioritization_out/`.

**Pasada en limpio (2026-09-04).** El árbol de §2 está creado: `comun/` con
`nucleo.py`, `tdr.py` y cinco tests (los tres instantáneos, en verde), los seis
notebooks de `createDB/` copiados más `05_qc_base.ipynb` y la documentación en
markdown, y los `funciones_<carpeta>.py` y notebooks de las tres carpetas de
análisis en el formato de §3. **Nada se ejecutó todavía**: `gon4/` está vacío y
ninguna corrida escribió resultados. El orden de §6 es el que hay que seguir para
llenarlo.

---

## 1. Cambios que le haría a tu enfoque

Nueve puntos, del más al menos importante. Los tres primeros cambian decisiones
de fondo; el resto son convenciones.

### 1.1 Un único lugar donde vive el modelo

Hoy la misma métrica está implementada tres veces con tres resultados distintos:
la pAUC del notebook (que hasta hoy tenía el signo invertido y sin
interpolación), `reproducir_v6/comun_tdr.py::auc01_mcclish`, y lo que cada
script de `reproducir_v6/` recalcula por su cuenta. El bug de McClish sobrevivió
meses justamente por eso.

**Propuesta.** En la raíz de `TDR_2026_v4/` una carpeta `comun/` con:

| Archivo | Qué contiene | Quién lo toca |
|---|---|---|
| `nucleo.py` | `get_druggable_targets`, `get_annot_druggability_pv`, `rs`, `nds` — el `nds_fun.py` actual, movido tal cual | **congelado**: se cambia sólo con un test que lo justifique |
| `tdr.py` | rutas, carga de tablas, métricas (`pauc`, `mcclish`, `auc01`), paralelización, paleta y estilo de figuras | compartido por las 4 carpetas |

Cada carpeta de análisis tiene además su propio `funciones_<carpeta>.py` con lo
que es específico de ese análisis. Regla: **si una función la necesitan dos
carpetas, sube a `comun/`; nadie reimplementa una métrica**.

*Por qué las carpetas no van numeradas:* un nombre que empieza con dígito no es
importable (`import 00_comun` es sintaxis inválida) y obligaría a trucos de
`sys.path` en cada notebook. El orden lo dan los números de los notebooks, y este
README.

### 1.2 El análisis corre en la celda, no en un script aparte

Con el barrido paralelizado (`n_core = 20`) el sweep completo son minutos, no
horas, y la desorfanización está en el mismo orden. A esa escala un punto de
entrada por línea de comandos sólo agrega indirección: conviene que **el ciclo
del análisis se lea en la celda de corrida**, que es donde uno lo va a mirar.

**Propuesta.** La celda 4 contiene el ciclo completo —sobre especies, sobre
drogas— y llama a las funciones de `funciones_<carpeta>.py` para cada paso. Las
funciones son las piezas; el orden en que se combinan se ve en el notebook:

```python
# celda de corrida: el ciclo se lee aca
for sp in spoi:
    ctx = fg.preparar_especie(datos, sp)                  # semilla LOSO + cat_rs
    ctx["rs"] = fg.relevance_por_alpha(ctx, params["alpha"])
    filas = fg.paralelizar(fg.evaluar_combinacion, fg.grilla_valida(params), n_core=n_core)
    fg.guardar_barrido(filas, SALIDAS, NB, params)        # -> 01_<sp>.csv
```

Lo que sí se mantiene del planteo anterior:

* **La celda de resultados no depende de la de corrida.** Abrir el notebook y
  correr 1, 2, 5, 6… produce todas las figuras leyendo de `gon4/`, sin
  re-ejecutar el análisis. Es lo que separa un notebook revisable de uno que hay
  que correr entero.
* **Nada de lógica dentro del ciclo.** Si la celda de corrida pasa de ~15 líneas,
  lo que sobra es una función que le falta al `.py`.
* Cada corrida deja su `meta.json` (§1.4) al lado de las salidas.

### 1.3 Se recalculan las priorizaciones: lo viejo usa la métrica equivocada

`genoma_completo_v3/` y `v4/` se calcularon con la pAUC invertida y sin
interpolación en el borde; `pseudohuerfanas/v3`, `v4` y `claude/` usan parámetros
óptimos leídos de esos CSV. Los números no son comparables entre sí ni con los
informes.

**Decidido.** Se **recalcula todo** dentro de v4: el barrido se vuelve a correr
desde `genome_prioritization/` hacia `gon4/genome_prioritization_out/`, y la
desorfanización se rehace con esos óptimos hacia `gon4/huerfanas_out/`.
Lo de `gon3/resultados/` queda como histórico y ninguna notebook de v4 lo lee;
`genoma_completo_v5/` sirve de control: los números nuevos tienen que coincidir
con los suyos.

### 1.4 Cada salida se acompaña de su `meta.json`

Hoy hay `v3`, `v4`, `v5`, `claude`, `it0`–`it3` y no hay forma de saber qué las
distingue sin leer el código que las produjo. Un `NN_meta.json` en cada
`<carpeta>_out/`, con el número del notebook que lo escribió: fecha, parámetros
de la corrida, fecha de las tablas de `DB/`, `n_core` y notebook de origen. Diez líneas de
código, y elimina la ambigüedad que hoy obliga a reconstruir la historia a mano.

### 1.5 `beta = 0` ya está en la grilla: no hace falta una carpeta aparte

Como la grilla es `beta ∈ {0, 0.5, 1.0, 1.5}`, la comparación **G′r (β=0) vs
G′rk (β>0)** ya está contenida en los 120 renglones de cada CSV del barrido. No
es un análisis nuevo: es una **vista** del barrido (`02_optimo_y_consistencia.ipynb`,
un panel de figura). Se elimina `red_gr_vs_grk/` como carpeta.

### 1.6 Sin comparación con la forma del RS del paper; α se queda

Queda fuera **sólo** `alpha_rs_paper/00_comparar_rs.py`, que contrasta la `rs()`
del código contra la forma funcional publicada. **α sigue siendo parámetro de la
grilla** y su efecto sobre la pAUC se analiza en `genome_prioritization/` igual
que el de β, λ y γ.

### 1.7 Migrar lo que entra al paper, y decirlo explícitamente

`reproducir_v6/` tiene 12 carpetas de scripts que en su mayoría **nunca se
corrieron completos**. Migrarlas todas a v4 traslada la deuda. La sección 4 de
este documento fija qué entra, qué entra recortado y qué se queda en v3.

### 1.8 Rutas en un solo lugar

`full_genome_v4.ipynb` apuntaba a `resultados/genoma_completo/`, un directorio
que ya no existe: la celda de visualización estaba rota y nadie se enteró. En v4,
todas las rutas son constantes de `comun/tdr.py` y las salidas se resuelven
contra el symlink `gon4`:

```python
V4   = Path("/home/ggiordano/TDR/TDR_2026_v4")
DB   = Path("/data1/TDR-2025/TDR_v7/gon3/DB")   # entrada; no se reescribe
GON4 = V4 / "gon4"                              # -> /data1/TDR-2025/TDR_v7/gon4

def out(carpeta):
    """Salidas del analisis `carpeta`: gon4/<carpeta>_out/ (la crea si no esta)."""
    d = GON4 / f"{carpeta}_out"
    (d / "figuras").mkdir(parents=True, exist_ok=True)
    return d
```

y en la celda 1 de cada notebook, `SALIDAS = tdr.out("genome_prioritization")`.
Cero literales de path dentro de los notebooks.

### 1.9 Numeración: en los notebooks y en sus salidas, no en las carpetas

Las carpetas quedan sin número, por lo del import. La numeración vive en los
notebooks (`01_`, `02_`, …) y **se propaga a todo lo que generan**:

```
genome_prioritization/01_barrido_parametros.ipynb
  -> gon4/genome_prioritization_out/01_sp21.csv … 01_sp4.csv
  -> gon4/genome_prioritization_out/01_meta.json
genome_prioritization/02_optimo_y_consistencia.ipynb
  -> gon4/genome_prioritization_out/02_optimos_por_especie.csv
  -> gon4/genome_prioritization_out/figuras/02_f01_grilla_beta.pdf (+ .png)
```

Así cada archivo dice de qué notebook salió, y las salidas de dos notebooks nunca
se pisan aunque tengan nombres parecidos.

Sobre el nombre en sí: **decidido** (2026-09-04) que quedan los actuales con el
typo corregido —`createDB`, `analiceDB`, `genome_prioritization`, `huerfanas`—,
que es lo que ya está escrito en las carpetas, en sus gemelas `_out` y en los
`sys.path`. La alternativa era pasar todo a castellano (`crear_db`,
`analizar_db`, `priorizacion_genoma`); se descartó por ser el cambio más caro
sobre lo que ya existe.

---

## 2. Estructura

```
TDR_2026_v4/
├── README.md                       ← este documento
├── gon4 -> /data1/TDR-2025/TDR_v7/gon4      todas las salidas van acá
│
├── comun/
│   ├── nucleo.py                   modelo (ex tdr-graph/nds_fun.py), congelado
│   ├── tdr.py                      rutas, carga, métricas, paralelización, estilo
│   └── tests/                      test_metricas, test_proyeccion, test_rutas,
│                                   test_integridad, test_fuga  (+ README)
│
├── createDB/                       copia de v3 (6 notebooks + documentación)
│   ├── 00_targets_species.ipynb … 04b_anotaciones_omcl.ipynb
│   ├── 05_qc_base.ipynb            NUEVO: chequeos de integridad de las 8 tablas
│   └── documentacion_createDB.md   el .txt actual, pasado a markdown
│
├── analiceDB/
│   ├── 01_dimensiones_red.ipynb
│   ├── 02_capas_y_conectividad.ipynb
│   ├── 03_calidad_bioactividades.ipynb
│   └── funciones_analice.py
│
├── genome_prioritization/
│   ├── 01_barrido_parametros.ipynb
│   ├── 02_optimo_y_consistencia.ipynb
│   └── funciones_genome.py
│
└── huerfanas/
    ├── 01_pseudohuerfanas.ipynb
    ├── 02_cobertura_semilla.ipynb
    ├── 03_desorfanizacion_global.ipynb
    ├── 04_aplicacion_sp26.ipynb
    └── funciones_huerfanas.py
```

Ninguna carpeta de análisis guarda resultados: cada una escribe en su gemela
`<carpeta>_out/` dentro de `gon4` (que es `/data1`, no el home), y todo archivo
lleva el número del notebook que lo generó.

```
gon4/
├── createDB_out/
│   ├── 05_qc_resumen.csv, 05_conteos_por_tabla.csv, 05_meta.json
│   └── figuras/  05_f01_qc.pdf|png
├── analiceDB_out/
│   ├── 01_dimensiones_capas.csv, 01_grado_afiliaciones.csv, 01_meta.json
│   ├── 02_*.csv, 03_*.csv
│   └── figuras/  01_f01_*.pdf|png, 02_f01_*, …
├── genome_prioritization_out/
│   ├── 01_<especie>.csv  (x16), 01_meta.json
│   ├── 02_optimos_por_especie.csv, 02_consistencia.csv
│   └── figuras/  01_f01_*, 02_f01_*, …
└── huerfanas_out/
    ├── 01_<especie>_pseudohuerfanas.csv, 01_meta.json
    ├── 02_cobertura_*.csv, 03_ranking_global.csv, 04_sp26_*.csv
    └── figuras/  01_f01_*, …
```

## 3. Formato de los notebooks

El que pediste, con nombres de celda fijos para que los cuatro análisis se lean
igual. Cada notebook es exactamente:

```
# <Título>                                          (markdown)

## Imports                                          (markdown)
[1] imports                                         una celda, sin lógica

## Datos                                            (markdown)
[2] carga                                           lee gon3/DB/ y salidas previas; sin transformar

## Acondicionamiento                                (markdown)
[3] preparación + vista de los datos de entrada     filtros, joins, y UNA figura/tabla
                                                    que muestra sobre qué se va a correr

## Corrida                                          (markdown)
[4] el ciclo del análisis                           llama a funciones_<X>.py; escribe NN_*.csv en <carpeta>_out/

# Resultados                                        (markdown)
[5] carga de lo generado                            lee de <carpeta>_out/, no recalcula
[6..n] figuras                                      una celda por figura, cada una guarda PDF+PNG
```

Reglas:

* **La celda 5 no depende de la 4.** Abrir el notebook y correr 1, 2, 5, 6… tiene
  que producir todas las figuras sin re-ejecutar el análisis. Es la diferencia
  entre un notebook que se puede revisar y uno que hay que correr entero.
* La celda 4 lleva **el ciclo a la vista** (§1.2): parámetros, `for` sobre
  especies o drogas, y una llamada por paso. Si pasa de ~15 líneas, lo que sobra
  es una función que le falta al `.py`.
* Cada celda de figura termina con `tdr.guardar(fig, f"{NB}_f01_<nombre>", FIGURAS)`
  y cada CSV se escribe como `SALIDAS / f"{NB}_<nombre>.csv"`: el número del
  notebook prefija todo lo que sale de él (§1.9).
* Nada se ejecuta al importar el `.py`.

### Plantilla de la celda 1 (idéntica en los cuatro)

```python
import sys, os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")          # antes de numpy: un hilo por proceso

import numpy as np, pandas as pd
import matplotlib.pyplot as plt, seaborn as sns

V4 = "/home/ggiordano/TDR/TDR_2026_v4"
sys.path.insert(0, f"{V4}/comun")
%load_ext autoreload
%autoreload 2
import tdr, nucleo as nf
import funciones_genome as fg              # el .py de esta carpeta

SALIDAS = tdr.out("genome_prioritization")   # -> gon4/genome_prioritization_out/
FIGURAS = SALIDAS / "figuras"
NB      = "01"                             # numero de este notebook: prefija todo lo que genere
n_core  = 20
```

---

## 4. Qué va en cada carpeta

### `createDB/` — copia de v3

Los 6 notebooks tal cual. Dos agregados baratos:

* `05_qc_base.ipynb`: los chequeos de `reproducir_v6/tests/test_integridad.py`
  vueltos notebook (conteos por tabla, huérfanos de FK, duplicados, targets sin
  especie), para poder afirmar que la base está sana antes de correr nada.
* Congelar en la documentación la **versión y fecha de ChEMBL/UniProt/InterPro**
  usadas. Hoy no está escrito en ningún lado.

### `analiceDB/` — descriptivos de la base y de la red

| Notebook | Contenido | De dónde sale |
|---|---|---|
| `01_dimensiones_red.ipynb` | nodos y enlaces de las 3 capas, disponibilidad por especie, distribución de grado de las afiliaciones | `analiceDB/01_tablas_analisis.ipynb` + `reproducir_v6/descriptivos_red/` |
| `02_capas_y_conectividad.ipynb` | tamaños de cluster, persistencia, perfiles hacia afuera, masas; conectividad entre especies por anotaciones y por drogas | `analiceDB/01` y `02` |
| `03_calidad_bioactividades.ipynb` | composición de las bioactividades por tag y fuente; filtro de promiscuidad química | `reproducir_v6/capa_quimica/02_filtro_promiscuidad.py` |

El **umbral de actividad (20 µM en el código vs 2 µM en el roadmap)** queda
**fuera del alcance de v4**: es una pregunta para Mercedes, no un análisis a
ejecutar acá. Se documenta el valor vigente y se deja anotado en §7.

Vale la pena mantener el punto 3: el filtro de promiscuidad ya mostró que **30
compuestos de MW < 150 Da concentran 104 075 relaciones de subestructura, el 12 %
de la capa**. Es un resultado que condiciona todo lo que viene después, y por eso
**en v4 el filtro se aplica** (§7.3): `03_calidad_bioactividades.ipynb` escribe
`03_compuestos_promiscuos.csv` y `huerfanas/` lo consume antes de armar ninguna
semilla. Con lo cual `analiceDB/03` deja de ser opcional: va **antes** del paso 5.

### `genome_prioritization/` — barrido y elección de parámetros

| Notebook | Contenido |
|---|---|
| `01_barrido_parametros.ipynb` | el `full_genome_v5.ipynb` actual: LOSO, 16 especies × 120 combinaciones, paralelizado con `n_core=20` (minutos). Salida `gon4/genome_prioritization_out/01_<especie>.csv`: un CSV por especie con `AUC01` (McClish + interpolación), `AUC`, `YoudenCutOff`, `N_targets` |
| `02_optimo_y_consistencia.ipynb` | ROC y distribución de scores de la especie elegida; efecto de cada parámetro sobre la pAUC (**α**, β, λ, γ); **β=0 vs β>0** (G′r vs G′rk); plateau del óptimo y estabilidad del top-K; consistencia entre especies y entre grupos taxonómicos; factor dominante |

Fuentes: `tdr-graph/full_genome_v5.ipynb`, `reproducir_v6/consistencia_parámetros/`,
`reproducir_v6/auc01_mcclish/` (ya absorbido en la métrica).

Opcional, si querés cerrar la Fase 5: un `03_priorizacion_multidimensional.ipynb`
con ortología humana + esencialidad + expresión (`reproducir_v6/priorizacion_multidimensional/`).
Lo dejaría para después de que el barrido y las huérfanas estén cerrados.

### `huerfanas/` — pseudohuérfanas, cobertura de semilla y desorfanización

| Notebook | Contenido | De dónde sale |
|---|---|---|
| `01_pseudohuerfanas.ipynb` | construcción de las k=1, grupos 0–3, `frank` por especie, clasificación de semilla (nula / no informativa / informativa) | `tdr-graph/orphan_drugs_v4.ipynb` |
| `02_cobertura_semilla.ipynb` | el cuello de botella: 61 % semilla nula. Umbral de similitud química, KEGG como tercera fuente, capa fenotípica de respaldo | `reproducir_v6/cobertura_semilla/` |
| `03_desorfanizacion_global.ipynb` | ranking global, `r*G`, distribución de `rSS`, partición directa/indirecta | `reproducir_v6/desorfanizacion_global/` + `inferencia_directa_indirecta/` |
| `04_aplicacion_sp26.ipynb` | embudo de compuestos huérfanos activos contra *la especie foco* y casos de estudio | `reproducir_v6/aplicacion_sp26/` |

`01` y `02` son el corazón: mientras el 61 % de las pseudohuérfanas tenga semilla
nula, ninguna reponderación mueve el número (las iteraciones it0–it3 dieron
resultados idénticos, que es lo que se espera si la limitación es estructural).

### Qué se queda en v3

| Carpeta | Motivo |
|---|---|
| `reproducir_v6/red_gr_vs_grk/` | subsumido por la grilla, que ya tiene β=0 (§1.5) |
| `reproducir_v6/alpha_rs_paper/` | decisión tuya: sin comparación de RS (§1.6) |
| `reproducir_v6/capa_quimica/00,01` | fingerprints 512/1024 y métricas alternativas: sólo si se decide regenerar la capa química |
| `reproducir_v6/umbral_bioactividad/` | el corte 20 µM vs 2 µM es una consulta a Mercedes, no un análisis de v4 (§7) |
| `gon3/resultados/genoma_completo_v3, v4` | métrica vieja; v5 se conserva como control (§1.3) |
| `gon3/resultados/pseudohuerfanas/claude/` | iteraciones it0–it3, sin efecto; el hallazgo ya está resumido en `02_cobertura_semilla` |
| presentaciones e informes de `reproducir_v6/` | son entregables, no análisis; se referencian desde acá |

---

## 5. Funciones: principales y auxiliares

El `.py` de cada carpeta se divide en dos bloques, con esta regla:

> **Principal** = cambia un número que va al paper (semilla, relevance score,
> propagación, métrica, ranking).
> **Auxiliar** = no lo cambia (I/O, paralelización, formateo, figuras, chequeos).

Cada principal lleva docstring con entradas/salidas **y un test** en
`comun/tests/`. Las auxiliares no necesitan test. El ciclo que las combina
**no** es una función: vive en la celda de corrida del notebook (§1.2).

### `genome_prioritization/funciones_genome.py`

```python
# ============================ PRINCIPALES ============================
def preparar_especie(datos, sp_code):
    """LOSO: semilla global sin `sp_code`, cat_rs, true positives y targets de la especie."""

def relevance_por_alpha(cat_rs, alphas):
    """RS para cada alpha de la grilla (una vez por especie, no por combinación)."""

def evaluar_combinacion(ctx, combo):
    """Una combinación (alpha, beta, lambda, gamma) -> AUC01, AUC, Youden, N_targets."""

def guardar_barrido(filas, salidas, nb, params):
    """Arma el DataFrame de resultados (alpha/beta/lambda/gamma legibles), lo ordena
    por AUC01 y lo escribe como <salidas>/<nb>_<especie>.csv + <nb>_meta.json."""

def optimos_por_especie(path_barrido, k=1):
    """Fila(s) óptima(s) por especie; k>1 devuelve el top-K para el análisis de plateau."""

# ============================= AUXILIARES ============================
def grilla_valida(params):        # (ia,ib,il,ig) sin las repeticiones de gamma
def paralelizar(func, items, n_core=20):   # ProcessPoolExecutor con fork
def cargar_barrido(path, spoi=None)        # concatena los CSV + columna especie
def escribir_meta(path, **campos)          # meta.json de la corrida
def fig_roc(...), fig_grilla_beta(...), fig_consistencia(...)
```

### `huerfanas/funciones_huerfanas.py`

```python
# ============================ PRINCIPALES ============================
def construir_pseudohuerfanas(datos)                  # k=1, grupos 0-3
def semilla_de_droga(datos, cid, ...)                 # CSN(m) pesado, sin fuga
def clasificar_semilla(seed, blanco, sta)             # nula / no informativa / informativa
def priorizar_droga(ctx, cid)                         # frank, rG, rSS, clase de inferencia
def guardar_huerfanas(filas, salidas, nb, params)     # <nb>_*.csv + <nb>_meta.json
def cobertura_por_palanca(datos, palanca, ...)        # umbral químico | KEGG | fenotipo

# ============================= AUXILIARES ============================
def muestra_pseudohuerfanas(k1, n=1000, rseed=123457)
def cargar_resultados(path), resumen_frank(df), escribir_meta(...)
def paralelizar(func, items, n_core=20)               # importada de tdr.py
def fig_grupos(...), fig_frank(...), fig_cobertura(...)
```

`analiceDB/funciones_analice.py` es casi todo auxiliar (conteos y figuras); su
único bloque principal son los descriptivos que se citan como cifras en el paper.

---

## 6. Orden de trabajo

1. **`comun/`** primero: mover `nds_fun.py` → `nucleo.py`, escribir `tdr.py`
   con rutas (`out()` hacia `gon4/<carpeta>_out/`), métricas y `paralelizar`, y
   correr los tests.
   Nada más se toca hasta que estén en verde.
2. **`createDB/`**: copia + `05_qc_base.ipynb`. Confirma que la base está sana.
3. **`genome_prioritization/`**: portar `full_genome_v5.ipynb` al formato de §3 y
   correr el barrido a `gon4/genome_prioritization_out/`. Control: los óptimos
   por especie tienen que coincidir con `gon3/resultados/genoma_completo_v5/`.
4. **`analiceDB/`**: independiente del barrido; se puede hacer en cualquier
   momento.
5. **`huerfanas/`**: `01` y `02` con los óptimos recién calculados; `03` y `04`
   después.

## 7. Decisiones que faltan

**Para resolver acá:**

1. ~~**Nombres de las carpetas**~~ — **decidido (2026-09-04)**: los actuales con
   el typo corregido (`genome_prioritization`). Ver §1.9.
2. **`priorizacion_multidimensional`**: ¿entra como `03_` dentro de
   `genome_prioritization/` o queda fuera del alcance de v4?
3. ~~**Filtro de promiscuidad química**~~ — **decidido (2026-09-04): se aplica.**
   Los 30 compuestos de MW < 150 Da que concentran 104 075 relaciones de
   subestructura (el **12 % de la capa**) se filtran antes de armar ninguna
   semilla. `analiceDB/03_calidad_bioactividades.ipynb` mide el efecto y escribe
   `03_compuestos_promiscuos.csv`; `huerfanas/` lo consume vía
   `tdr.compuestos_promiscuos()` y `tdr.filtrar_capa_quimica()`. **Consecuencia:**
   los números de `huerfanas/` no son comparables con los de
   `gon3/resultados/pseudohuerfanas/`, y `analiceDB/03` tiene que correr antes
   que el paso 5 del orden de trabajo.
4. **`gon4` y git**: los resultados pesados viven en `/data1` vía symlink. Falta
   decidir si `TDR_2026_v4/` va a git y con qué `.gitignore` (v3 no es un repo).
   El symlink `gon4` nunca debe entrar al repo.

**Para preguntarle a Mercedes** (no se ejecuta en v4):

5. **Umbral de bioactividad**: el código usa **20 µM** para etiquetar `positive`;
   el roadmap menciona **2 µM**, y el paper 2019 (Table 4) reporta con otro
   criterio. Hay que confirmar cuál es el corte correcto y si se aplica por tipo
   de ensayo (IC50/Ki/EC50) antes de volver a etiquetar nada. Si el corte cambia,
   se rehace `03_bioactivities_target_compound.csv` y con eso **todo** el árbol:
   barrido, óptimos y huérfanas. Por eso queda como consulta previa y no como
   análisis pendiente.
