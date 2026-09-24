# Plan de trabajo — proyecto final GW-AI

> **Actualización del 2026-09-24:** se portó desde `TDR_2026_v4` la corrección
> de un error del modelo en la semilla de desorfanización y se rehízo
> `huerfanas/` (registro en `historia/README_GWAI_correccion_semilla.md`).
>
> **Estado al 2026-09-10.** Las seis fases están hechas. Queda pendiente la
> publicación de la página (se espera el contacto de los organizadores) y un
> puñado de mejoras que se listan al final. El plan original, del 2026-09-08,
> está en el historial de git (`git show 0ab3fc3:PLAN.md`); acá cada ítem dice
> qué se pensaba hacer, qué se hizo y dónde está la evidencia.

Marcas: ✅ hecho · 🔀 hecho de otra forma que la planeada (se explica por qué)
· ⏳ pendiente.

---

## Parte A — Portar los análisis de v4 a la base codificada

La base nueva es toda enteros. Eso no rompe con una excepción: rompe en
silencio (un filtro contra el texto `"positive"` devuelve cero filas y el
análisis sigue). Por eso las pruebas iban antes que las corridas.

| # | qué se pensaba hacer | estado | evidencia |
|---|---|---|---|
| A.0 | Copiar de v4 `comun/`, los tres análisis y el informe, sin `createDB/` ni bibliografía | ✅ | commit `fe2b836`. El informe de v4 (`paper/`) después salió del repo: nombraba las especies (ver decisiones, D.4) |
| A.1 | Que el filtro de positivos no quede vacío con la base codificada | ✅ | `TAG_POSITIVE = 2` en `comun/tdr.py`; `test_el_filtro_de_positivos_no_quedo_vacio` |
| A.2 | Traducir los literales de texto (tags, tipo InterPro, códigos de especie, faltante, pseudo-especie) | ✅ | constantes comentadas en `comun/tdr.py`; PROVENANCE §3.1 |
| A.3 | Confirmar que sobreviven las 16 especies del modelo | ✅ | `test_las_16_especies_estan_completas` |
| A.4 | Rutas relativas a la raíz del repo; control de v5 copiado | ✅ | `tdr.RAIZ`; `test_todo_cuelga_del_repositorio`; `control/genoma_completo_v5/` |
| A.4 | Publicar los diccionarios chicos en `DB/mapeos/` + `comun/codigos.py` | 🔀 | no se publicó ningún diccionario: sólo las cinco constantes que el modelo necesita. Publicarlos habría deshecho la codificación (PROVENANCE §3.1) |
| A.4 | Resolver `RAW` (datos crudos) | 🔀 | ningún notebook los lee: el único cálculo que los necesita, la lista de compuestos promiscuos, pasó a `datos_externos/promiscuidad/derivar.py` y su salida se versiona (PROVENANCE §3.13) |
| A.5 | Leer la capa de aristas partida en dos `.gz` y desescalar el peso | ✅ | `tdr.leer_aristas()`: 36 152 622 aristas |
| A.6 | `cluster_consistent` como booleano | ✅ | cast explícito en `cargar_db` |
| A.7 | Explotar la asimetría: lo que no toca la capa química no debe moverse | ✅ | `verificacion/10`: priorización idéntica, 0 discrepancias. Con una corrección: `huerfanas/` sí toca la capa química y sí se mueve (PROVENANCE §5) |
| A.8 | Que InterPro y OrthoMCL no se pisen (hallazgo de las pruebas) | ✅ | prefijo repuesto en la carga; `test_los_dos_vocabularios_no_se_pisan` |
| A.9 | Cambiar una sola línea de `nucleo.py` (clave del diccionario de salida) | ✅ | PROVENANCE §3.4 |

---

## Parte B — Completar el proyecto según la consigna

La consigna (`final-project.html`) pide procedencia de cada resultado, una
verificación detrás de cada figura, una página HTML y un PDF, y un repositorio
que un agente pueda reproducir sin haber hablado con nadie.

### Fase 1 — Reproducibilidad de la infraestructura ✅

| qué se pensaba hacer | estado | evidencia |
|---|---|---|
| Raíz derivada del archivo, cero rutas absolutas | ✅ | `tdr.RAIZ`; CONVENCIONES §6 |
| `cargar_db()` adaptado a la base codificada | ✅ | Parte A |
| `requirements.txt` con versiones exactas | ✅ | entorno `TDR_GW`, python 3.12.11, pandas 3.0.5; `test_requirements_cubre_lo_que_se_usa` |
| Copiar el control de v5 | ✅ | `control/genoma_completo_v5/` |

### Fase 2 — Pruebas que atrapen la falla silenciosa ✅

| qué se pensaba hacer | estado | evidencia |
|---|---|---|
| Positivos no vacíos, 16 especies, aristas completas, base entera, métricas | ✅ | `comun/tests/`, 60 pruebas; la primera corrida encontró 6 fallas reales |
| (agregado) Que el entorno tenga todo lo que se importa | ✅ | `test_entorno.py` |
| (agregado) Que los resultados reproduzcan el control | ✅ | `test_resultados.py` |
| (agregado) Que las figuras se dibujen sin tocar la base | ✅ | `test_figuras.py` |

### Fase 3 — Correr los notebooks y comparar contra el oráculo ✅

| qué se pensaba hacer | estado | evidencia |
|---|---|---|
| Correr los 10 notebooks en orden | ✅ | corridas del 2026-09-08 y 2026-09-10; tablas y `NN_meta.json` en `resultados/` |
| Notebook de equivalencia contra v4 que clasifique cada diferencia | ✅ | `verificacion/10_equivalencia.ipynb`: 353 columnas, 252 idénticas, 97 explicadas por el recorte, 0 discrepancias |
| (agregado) Análisis de las pseudohuérfanas con semilla informativa | ✅ | `huerfanas/03`, figura `03_f05` |

### Fase 4 — Procedencia explícita ✅

| qué se pensaba hacer | estado | evidencia |
|---|---|---|
| Por cada resultado: qué lo produjo, qué entró, qué es propio y qué librería, qué alternativa se descartó | ✅ | `PROVENANCE.md` |
| Por cada figura: notebook, tablas, verificación | ✅ | `FIGURAS.md`, generado del registro de `comun/figuras.py` |
| `meta.json` por corrida | ✅ | `resultados/*/NN_meta.json` |
| Traer al repo las decisiones de `acondicionarDB/` | 🔀 | resumidas en PROVENANCE §2 y en `DB/meta.json`; `acondicionarDB/` no entra al repo porque lee la base real |
| (agregado) Separar la corrida de las figuras | ✅ | `comun/figuras.py`; CONVENCIONES §2 |

### Fase 5 — Los entregables para personas ✅ (publicación ⏳)

| qué se pensaba hacer | estado | evidencia |
|---|---|---|
| PDF | 🔀 | no se actualizó el informe de v4: se escribió uno nuevo, `informe/informe.pdf`, con todas las cifras leídas de las tablas (D.3) |
| Página HTML | ✅ | `index.html`, generada por `informe/generar_pagina.py` |
| (agregado) Presentación de 5 minutos | ✅ | `informe/presentacion.pdf` |
| Publicar la página y enlazarla desde la página general del curso | ⏳ | se espera el contacto de los organizadores |

### Fase 6 — Prueba de fuego ✅

| qué se pensaba hacer | estado | evidencia |
|---|---|---|
| Clonar el repo en un directorio limpio y reproducir siguiendo sólo el README | ✅ | hecho el 2026-09-10. Encontró que sin las tablas no se podía verificar nada y que la receta de entregables los pisaba con «??»: se versionaron las tablas y se agregaron guardias. `make verificar` pasa en un clon limpio |
| Corrida completa desde un clon | ⏳ | no repetida desde un clon (~45 min); sí completa en el repositorio de trabajo el 2026-09-10 |

---

## Decisiones que estaban abiertas

| # | pregunta | decisión |
|---|---|---|
| D.1 | `nombres_ipr()`: sin el diccionario privado no hay nombres de dominio | ✅ las familias se presentan por código; la función falla explicando por qué (PROVENANCE §3.6). Con una excepción posterior: los 391 dominios de quinasa sí se nombran, porque `genome_prioritization/03` los necesita (PROVENANCE §3.15) |
| D.2 | Alcance: ¿los tres análisis o uno a fondo? | ✅ los tres, con la priorización y la desorfanización como resultados principales |
| D.3 | ¿Se actualiza el informe de v4 o se escribe uno nuevo? | ✅ uno nuevo, más corto, con el contraste con v4 en una sola sección |
| D.4 | `paper/` nombraba las especies | ✅ salió del repo (commit `52a09dc`); sus tablas quedaron, anonimizadas, en `oraculo_v4/` |

---

## Lo que queda

| qué | por qué importa |
|---|---|
| ⏳ Publicar la página | la consigna la enlaza desde la página del curso |
| ⏳ Un control independiente del **modelo**, no sólo del port | el control de v5 corre el mismo `nucleo.py` (CONVENCIONES §3) |
| ✅ Que ningún notebook necesite datos privados | `analiceDB/03` ya no lee datos crudos; la lista de promiscuos es un insumo externo versionado (PROVENANCE §3.13) |
| ✅ Revisar por qué la semilla nula bajó del 83 al 69 % respecto de v4 | era el error de la semilla (vecinos buscados con el id del compuesto en tablas de clusters), que la codificación hacía caer en otros clusters. Corregido el 2026-09-24 con `comun/tests/test_vecinos.py`; con la corrección, las semillas coinciden con las de v4 corregida (PROVENANCE §3.16 y §5) |
| ✅ La lista de promiscuos tenía filas duplicadas (7 filas, 5 compuestos; en v4, 30 y 25) | `derivar.py` deduplica antes del cruce; el filtro ya la usaba como conjunto, así que ninguna semilla cambió (`datos_externos/promiscuidad/README.md`) |
| ⏳ La sensibilidad de r*G, más allá de k_σ | `03_rg_sensibilidad.csv` barre k_σ entre 2 y 4 (68–88); en v4 se vio que con toda la grilla (ventana, `desde`, suavizado) va de 20 a 195. Se presenta como corte operativo |
| ⏳ Palancas reales para ampliar la semilla | bajar el umbral de similitud por debajo de 0.8 exige recalcular la capa química; KEGG exige un mapeo que no está |
| ⏳ Verificar contra el artículo la partición directa/indirecta 68/32 | hoy se cita desde un comentario del código de v3 |
