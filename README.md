# TDR-GW-AI-course

Proyecto final del curso GW-AI: priorización de blancos terapéuticos sobre una
red multicapa, reproducida sobre una base codificada y con la procedencia de
cada resultado registrada.

> **In English.** Given a genome, rank its proteins by how promising they are as
> drug targets, by propagating scarce experimental evidence across a three-layer
> network (compounds, proteins, functional annotations). Evaluated leaving out
> all evidence of the organism being ranked, it beats chance in all 16 organisms;
> the same model recovers the target of single-target compounds whenever their
> chemical neighbourhood shares an annotation with it. Everything runs on an
> integer-coded version of the database, every number in the deliverables is read
> from the result tables, and every figure is redrawn from versioned tables in
> ~40 s (`make figuras`). The report and page are in Spanish; `CORRECCION.md`
> maps each requirement of the assignment to its evidence.

**Para corregir, empezar por [`CORRECCION.md`](CORRECCION.md).**

## Entregables

| archivo | qué es |
|---|---|
| `index.html` | la página de presentación |
| `informe/informe.pdf` | el informe (8 páginas) |
| `informe/presentacion.pdf` | la presentación de 5 minutos |

Los tres se arman leyendo `resultados/`: ninguno tiene una cifra escrita a mano.

## Instalar

```bash
conda create -n TDR_GW python=3.12
conda activate TDR_GW
pip install -r requirements.txt
```

Compilar el informe y la presentación pide además `pdflatex` y `bibtex`.

## Verificar sin correr nada (~1 min)

Las tablas de todas las corridas están **versionadas** en
`resultados/<analisis>_out/`, junto con el `NN_meta.json` de cada una. Con
ellas, **todas las figuras se reproducen desde las tablas pregeneradas**, sin
cargar la base ni correr ningún modelo:

```bash
make verificar      # 63 pruebas + las 36 figuras desde las tablas + cifras del informe
make figuras        # sólo las figuras, ~40 s
make entregables    # figuras, FIGURAS.md, informe, presentación y página, ~1,5 min
```

`make verificar` corre las pruebas, redibuja cada figura y comprueba que las
cifras que `informe/numeros.py` lee hoy de las tablas sean las mismas que están
en el informe versionado. Qué tablas lee cada figura y qué la verifica está en
[`FIGURAS.md`](FIGURAS.md).

Una figura sola, o las de un análisis:

```bash
python comun/figuras.py huerfanas
python comun/figuras.py 02_f01_roc_26
python comun/figuras.py --lista           # figura -> tablas que lee
```

## Tiempo de cómputo

Reproducir las tablas es lo único caro. Medido en la corrida del 2026-09-10, en
una máquina de 48 núcleos con **20 procesos**:

| notebook | qué hace | tiempo |
|---|---|---|
| `analiceDB/01_dimensiones_red` | dimensiones de la red, disponibilidad por especie | 21 s |
| `analiceDB/02_capas_y_conectividad` | tamaños de cluster, componentes conexas, conectividad | 3,5 min |
| `analiceDB/03_calidad_bioactividades` | bioactividades por tag, efecto del filtro de promiscuidad sobre la base | 1 min |
| `genome_prioritization/01_barrido_parametros` | 16 especies × 120 combinaciones | 7,5 min |
| `genome_prioritization/02_optimo_y_consistencia` | óptimo, meseta, consistencia | 33 s |
| `genome_prioritization/03_comparacion_beta_promiscuidad` | β = 1 contra β = 0 a parámetros fijos; dónde quedan las quinasas | 2,5 min |
| `huerfanas/01_pseudohuerfanas` | 7 782 compuestos, uno por uno | 13 min |
| `huerfanas/02_cobertura_semilla` | tres palancas sobre 500 compuestos | 17 min |
| `huerfanas/03_desorfanizacion_global` | r*G, directa/indirecta, semilla informativa | 21 s |
| `huerfanas/04_aplicacion_especie` | embudo de la especie foco | 1 min |
| `verificacion/10_equivalencia` | comparación contra `oraculo_v4/` | 10 s |
| **total** | | **~48 min** |

La carga más pesada de la base (con la capa química, 36 M aristas) tarda ~30 s
y ocupa 2,1 GB; los procesos de la corrida paralela la comparten por `fork`. La
cantidad de procesos se cambia con `TDR_NCORE` (por defecto 20); con menos, los
tres notebooks largos tardan proporcionalmente más.

En cambio, **redibujar las 36 figuras desde las tablas lleva ~40 s**, y las
pruebas, ~25 s.

## Reproducir la corrida (~45 min)

```bash
make corrida        # los 10 notebooks en orden (CONVENCIONES.md §9)
make entregables
```

Cada notebook escribe sus tablas y su `NN_meta.json` en
`resultados/<analisis>_out/`. Se puede correr uno solo con
`jupyter nbconvert --to notebook --execute --inplace <notebook>` desde su
carpeta. Todos tienen las mismas secciones: una **corrida** que carga la base y
escribe tablas, y **Resultados**, que sólo las lee y dibuja; con la celda de
imports y esa sección alcanza para ver las figuras sin volver a correr nada
(CONVENCIONES.md §2).

### El único insumo calculado fuera: la lista de compuestos promiscuos

Todos los notebooks corren desde `DB/`. La excepción es la **lista** de
compuestos promiscuos que filtra `huerfanas/` (CONVENCIONES.md §10): contar
cuántas superestructuras contienen a cada compuesto exige las relaciones de
subestructura crudas, que `DB/` no conserva. Esa lista se calcula aparte, con
`datos_externos/promiscuidad/derivar.py`, sobre tres archivos privados (datos
crudos de ChEMBL y el diccionario de la codificación), y se versiona como
insumo, igual que `DB/`. `analiceDB/03` la lee y mide desde la base cuánto saca.
Qué archivos necesita, cómo se produjo y cómo se verificó está en
[`datos_externos/promiscuidad/README.md`](datos_externos/promiscuidad/README.md).

## Datos

`DB/` es la base acondicionada: 260 MB, sólo enteros. La generó
`acondicionarDB/acondicionar.py` (fuera del repositorio) a partir de la base
original, en tres operaciones: **recorte** de las componentes químicas sin
bioactividad positiva, **codificación** a enteros con semilla fija y
**empaquetado** de la capa de aristas en dos `.csv.gz`. El registro completo
(semilla, criterio de recorte, diccionarios, esquema de cada tabla) está en
`DB/meta.json`, y las decisiones con su alternativa descartada, en
`PROVENANCE.md` §2.

Las equivalencias que el modelo necesita están escritas en `comun/tdr.py`
(`TAG_POSITIVE`, `IPR_DOMAIN`, `SPECIES`, …), con el comentario de qué
selecciona cada una. Las especies **se nombran**: `SPECIES` da, para cada uno de
los 16 códigos, el organismo (`26` → *Plasmodium falciparum*) junto con su tipo
(grupo, reino, si es parásito), y las figuras y las tablas de `resultados/` se
rotulan con el binomio abreviado. El código sigue siendo la clave con la que la
especie entra en `DB/`, en `control/` y en las tablas.

El resto de la codificación sigue en pie: compuestos, blancos, clusters,
InterPro y OrthoMCL no se traducen, y los diccionarios de `mapeos/` no se
publican. Eso es lo que impide reidentificar los compuestos (`PROVENANCE.md`
§3.13); los nombres de las 16 especies no aportan nada a esa reidentificación.

## Estructura

| carpeta o archivo | qué es |
|---|---|
| `CORRECCION.md` | la consigna, requisito por requisito, y dónde está la evidencia |
| `FIGURAS.md` | cada figura: notebook, tablas que lee, dónde aparece y qué la verifica |
| `PROVENANCE.md` | de dónde viene cada resultado y qué alternativa se descartó |
| `CONVENCIONES.md` | las reglas de organización que cita el código |
| `PLAN.md` | el plan de trabajo, con lo que se hizo y lo que no |
| `DB/` | la base codificada |
| `comun/` | `tdr.py` (rutas, carga, métricas, estilo), `nucleo.py` (el modelo), `figuras.py` (registro de figuras), `tests/` |
| `analiceDB/`, `genome_prioritization/`, `huerfanas/` | los tres análisis, 10 notebooks; cada carpeta con `funciones_*.py` (cálculo) y `figuras_*.py` (figuras) |
| `verificacion/` | la comparación contra el oráculo |
| `resultados/` | las tablas de la corrida, versionadas; las figuras se regeneran |
| `control/` | óptimos por especie de una corrida anterior, para contrastar |
| `oraculo_v4/` | las mismas tablas calculadas antes del recorte y de la codificación |
| `datos_externos/` | lo que no sale de `DB/`: la lista de compuestos promiscuos (insumo versionado) y la procedencia de `tdr.ANN_QUINASA` |
| `informe/` | fuentes del informe y de la presentación, y el generador de la página |
| `historia/` | el documento de trabajo de la versión v4, sólo como registro |
| `final-project.html` | la consigna del curso |
