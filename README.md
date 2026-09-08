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
python -m unittest discover -s comun/tests -p "test_*.py"    # 38 tests, ~15 s
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
| `PROVENANCE.md` | de dónde viene cada resultado y qué alternativa se descartó |
| `DB/` | la base codificada; sólo enteros |
| `comun/` | `tdr.py` (rutas, carga, métricas, figuras), `nucleo.py` (el modelo), `tests/` |
| `analiceDB/`, `genome_prioritization/`, `huerfanas/` | los tres análisis, 9 notebooks |
| `verificacion/` | la comparación contra el oráculo |
| `resultados/` | salidas de las corridas (el contenido no se versiona) |
| `control/` | óptimos por especie de una corrida independiente, para contrastar |
| `oraculo_v4/` | las mismas tablas calculadas antes del recorte y de la codificación |

## Regenerar los entregables

```bash
python informe/numeros.py          # los números del informe, desde las tablas
make -C informe                    # -> informe/informe.pdf
python informe/generar_pagina.py   # -> index.html
```

Ninguno de los dos entregables tiene un número escrito a mano: los dos se
arman leyendo `resultados/`.
