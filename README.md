# TDR-GW-AI-course

Proyecto final del curso GW-AI.

## Contenido

- `final-project.html` — consigna y página de presentación del proyecto.
- `DB/` — base de datos acondicionada del proyecto (CSV y CSV comprimidos).
  - `DB/meta.json` describe cada tabla: filas, columnas, tamaño, y el criterio
    de recorte aplicado al generarla.

## Datos

`DB/` fue generado por `acondicionar.py` a partir de la DB original de TDR.
El detalle completo (semilla, escala de peso, criterio de recorte, diccionarios
y esquema de cada tabla) está en `DB/meta.json`.

`01_edges_clusters_fingerprint` está partido en dos archivos `.gz`
(`part00`, `part01`) que deben concatenarse para reconstruir la tabla completa.

## Cómo correr los análisis

```bash
conda create -n TDR_GW python=3.12
conda activate TDR_GW
pip install -r requirements.txt
python -m unittest discover -s comun/tests -p "test_*.py"    # 51 tests, ~25 s
```

Los notebooks se corren en este orden, y cada uno escribe en
`resultados/<carpeta>_out/` junto con su `NN_meta.json` de procedencia:

| orden | carpeta | qué hace |
|---|---|---|
| 1 | `analiceDB/` | descriptivos de la base y de la red (3 notebooks) |
| 2 | `genome_prioritization/` | barrido de parámetros y elección del óptimo (2) |
| 3 | `huerfanas/` | pseudohuérfanas, cobertura de semilla, desorfanización (4) |

La celda 1 de cada notebook agrega `comun/` al path; no hace falta configurar
nada más, porque todas las rutas se derivan de la ubicación del repositorio.

### Corrida y figuras, por separado

Cada notebook tiene dos mitades que cuestan muy distinto:

| mitad | qué hace | cuánto tarda |
|---|---|---|
| **corrida** (Datos → Acondicionamiento → Corrida) | carga la base, calcula y escribe las tablas en `resultados/<analisis>_out/` | de segundos a ~15 min por notebook |
| **Resultados** | sólo lee esas tablas y dibuja con `F.dibujar(...)` | segundos |

Las figuras no se dibujan en ninguna celda de corrida: están todas en
`<carpeta>/figuras_<analisis>.py`, cada una declarando qué tablas lee, y un
test (`comun/tests/test_figuras.py`) verifica que ninguna toque la base. Para
retocar o regenerar figuras no hace falta volver a correr nada:

```bash
python comun/figuras.py                   # todas las figuras, ~30 s
python comun/figuras.py huerfanas         # las de un análisis
python comun/figuras.py 02_f01_roc_26     # una sola
python comun/figuras.py --lista           # qué tablas lee cada figura
```

Desde un notebook alcanza con correr la celda de imports y saltar a la sección
«Resultados».

## Sobre los códigos de la base

`DB/` está codificada: toda columna es un entero y los diccionarios que traducen
esos enteros a la notación original son privados. Las equivalencias que el
modelo necesita están escritas en `comun/tdr.py` (`TAG_POSITIVE`, `IPR_DOMAIN`,
`SPECIES`, …), con el comentario de qué selecciona cada una.

Las especies se identifican por código, no por nombre. **El nombre de la especie
no aparece en ninguna parte del repositorio**: ni en el código, ni en los datos,
ni en los nombres de archivo, ni en la documentación. Lo que sí se publica es el
tipo de organismo (grupo, reino, si es parásito), porque el modelo lo usa y sin
él las figuras por grupo no se leen. Tres de las 16 son la única de su
combinación (grupo, parásito), lo que para ellas equivale a identificarlas: está
anotado en `comun/tdr.py`.

## Estructura

| carpeta | qué es |
|---|---|
| `index.html` | **la página de presentación**, generada desde las tablas |
| `informe/informe.pdf` | **el informe**, con sus números generados desde las tablas |
| `informe/presentacion.pdf` | **la presentación** de 5 minutos, con las mismas cifras y figuras |
| `PROVENANCE.md` | de dónde viene cada resultado y qué alternativa se descartó |
| `DB/` | la base codificada; sólo enteros |
| `comun/` | `tdr.py` (rutas, carga, métricas, estilo), `nucleo.py` (el modelo), `figuras.py` (registro y regeneración de figuras), `tests/` |
| `analiceDB/`, `genome_prioritization/`, `huerfanas/` | los tres análisis, 9 notebooks; cada carpeta con `funciones_*.py` (cálculo) y `figuras_*.py` (figuras) |
| `verificacion/` | la comparación contra el oráculo |
| `resultados/` | salidas de las corridas (el contenido no se versiona) |
| `control/` | óptimos por especie de una corrida independiente, para contrastar |
| `oraculo_v4/` | las mismas tablas calculadas antes del recorte y de la codificación |

## Regenerar los entregables

```bash
python informe/numeros.py          # los números del informe, desde las tablas
make -C informe                    # redibuja las figuras y compila -> informe/informe.pdf
                                   #   y la presentacion de 5 min -> informe/presentacion.pdf
python informe/generar_pagina.py   # -> index.html
```

Ninguno de los dos entregables tiene un número escrito a mano: los dos se
arman leyendo `resultados/`.
