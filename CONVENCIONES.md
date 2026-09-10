# Convenciones del código

Las reglas que sostienen la organización del repositorio. El código las cita
como `CONVENCIONES.md §N`. Vienen de la reorganización de la versión v4 del
proyecto (su documento de trabajo está en `historia/README_TDR_2026_v4.md`, sólo
como registro), y acá están escritas tal como rigen en este repositorio.

## §1 · Un único lugar para el modelo

`comun/nucleo.py` es el modelo (semilla, relevance score, propagación) y está
**congelado**: se cambia sólo con una prueba de `comun/tests/` que lo justifique.
En este repositorio se cambió una línea, documentada en `PROVENANCE.md` §3.4.
`comun/tdr.py` tiene lo compartido: rutas, carga, métricas, paralelización y
estilo. Cada carpeta de análisis tiene su `funciones_<analisis>.py` (cálculo) y
su `figuras_<analisis>.py` (figuras). **Si una función la necesitan dos carpetas,
sube a `comun/`; nadie reimplementa una métrica.** El motivo: en la versión v3
la misma métrica estaba implementada tres veces con tres resultados distintos, y
un error en una de ellas sobrevivió meses.

## §2 · Formato de los notebooks

Todos tienen las mismas secciones, en este orden:

| sección | qué hace |
|---|---|
| Imports | una celda, sin lógica |
| Datos | carga la base y las salidas previas, sin transformar |
| Acondicionamiento | filtros, joins y una vista de los datos de entrada |
| Corrida | **el ciclo del análisis a la vista**: llama a `funciones_*.py` y escribe tablas `NN_*.csv` en `resultados/<analisis>_out/` |
| Resultados | sólo lee esas tablas y dibuja con `F.dibujar(...)` |

La sección Resultados **no depende de la corrida**: con la celda de imports
alcanza. Ninguna celda dibuja por su cuenta; las figuras viven en
`figuras_*.py` y leen sólo tablas declaradas (`comun/figuras.py`). Lo vigila
`comun/tests/test_figuras.py`.

## §3 · Control externo

Los óptimos por especie se contrastan contra `control/genoma_completo_v5/`, una
corrida anterior hecha sobre la base sin codificar y con la métrica ya corregida.
Tienen que coincidir exactamente (`comun/tests/test_resultados.py`). El control
corre el mismo modelo: valida el port, no el modelo.

## §4 · Procedencia de cada corrida

Toda celda de corrida termina escribiendo `NN_meta.json` junto a sus tablas
(`tdr.escribir_meta`): fecha, notebook, parámetros, fecha de cada tabla de `DB/`,
número de procesos y versiones de python y pandas.

## §5 · β = 0 está en la grilla

La grilla del barrido incluye `beta ∈ {0, 0.5, 1.0, 1.5}`, así que la
comparación G′r (β = 0) contra G′rk (β > 0) es una vista del barrido
(`genome_prioritization/02`), no un análisis aparte.

## §6 · Rutas

Todo se resuelve contra la raíz del repositorio (`tdr.RAIZ`, la carpeta que
contiene `DB/`). Ningún notebook tiene un literal de ruta. `DB/` es de sólo
lectura; las salidas van a `resultados/<analisis>_out/`, que crea `tdr.out()`.
Ningún notebook lee nada fuera del repositorio. El único cálculo que necesita
datos privados, la lista de compuestos promiscuos, vive en
`datos_externos/promiscuidad/derivar.py` (rutas `TDR_RAW` y `TDR_MAPEOS`) y su
salida se versiona como insumo (§10).

## §7 · Numeración

Las carpetas no llevan número, así son importables desde Python. Los notebooks
sí, y **todo archivo que un notebook genera arranca con su número**: el `01`
escribe `01_*.csv`, `01_meta.json` y `01_fNN_*.pdf`. Mirando una tabla o una
figura se sabe qué notebook la produjo.

## §8 · Funciones principales y auxiliares

Cada `funciones_*.py` separa dos bloques. **Principal**: cambia un número que va
al informe (semilla, relevance score, propagación, métrica, ranking); lleva
docstring con entradas y salidas. **Auxiliar**: no lo cambia (lectura y
escritura, paralelización, formateo). El ciclo que las combina no es una
función: vive en la celda de corrida (§2).

## §9 · Orden de corrida

```
analiceDB/01, 02, 03  ─┐
genome_prioritization/01 → 02 ─┤
                                └─> huerfanas/01 → 02, 03 → 04  ─>  verificacion/10
```

`huerfanas/` necesita los óptimos de `genome_prioritization/02` y la lista de
compuestos promiscuos de `analiceDB/03` (§10). `verificacion/10` compara todo
contra `oraculo_v4/`.

## §10 · El filtro de promiscuidad química se aplica

Los compuestos chicos y promiscuos (peso molecular < 150 Da y más de 100
superestructuras que los contienen, el criterio del trabajo de 2016) se retiran
de la capa de subestructuras **antes** de armar cualquier semilla.
La lista sale de `datos_externos/promiscuidad/` (se calcula desde datos crudos,
fuera de los análisis); `analiceDB/03` la transcribe a
`03_compuestos_promiscuos.csv` y mide desde la base cuánto saca, y `huerfanas/`
la aplica con `tdr.filtrar_capa_quimica()`.
