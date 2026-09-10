# Índice de figuras

Generado por `python comun/figuras.py --indice`; no editar a mano.

Cada figura se dibuja **sólo** desde las tablas de `resultados/<analisis>_out/`
que figuran en la columna *lee*, que están versionadas: `python comun/figuras.py`
las regenera todas en unos 40 segundos sin cargar la base. Cada tabla lleva el
número del notebook cuya celda de corrida la escribió (CONVENCIONES.md §7), y
el `NN_meta.json` de ese notebook registra la corrida.

La columna *verificación* dice qué respalda cada figura. Cuando no hay control
independiente, lo dice.

## analiceDB

Funciones en `analiceDB/figuras_analice.py`.

| figura | notebook | lee | aparece en | verificación |
|---|---|---|---|---|
| `01_f01_disponibilidad`<br>disponibilidad() | `01_dimensiones_red.ipynb` | `01_disponibilidad_especie.csv` | informe, página | Tabla idéntica a la del oráculo v4 (`verificacion/10`); `test_resultados` fija el total de druggables. |
| `01_f02_grado_afiliaciones`<br>grado() | `01_dimensiones_red.ipynb` | `01_grado_afiliaciones.csv` | — | `test_resultados.test_la_capa_de_anotaciones_no_se_movio`: categorías y afiliaciones idénticas a las de antes del recorte. |
| `01_f03_tramos_de_grado`<br>tramos_de_grado() | `01_dimensiones_red.ipynb` | `01_grado_afiliaciones.csv` | — | La misma tabla que `01_f02`, con el mismo control. |
| `02_f01_tamanos_cluster`<br>tamanos_cluster() | `02_capas_y_conectividad.ipynb` | `02_tamanos_cluster.csv` | — | Sin control externo: describe la capa química tal como quedó tras el recorte, cuyos conteos registra `DB/meta.json`. |
| `02_f02_componentes_conexas`<br>componentes_conexas() | `02_capas_y_conectividad.ipynb` | `02_componentes_conexas.csv` | — | Es una verificación del recorte en sí: las dos series tienen que coincidir, porque no debe quedar ninguna componente sin bioactividad positiva. |
| `02_f03_conectividad_anotaciones`<br>conectividad_anotaciones() | `02_capas_y_conectividad.ipynb` | `02_conectividad_anotaciones.csv` | — | Sin control: contra el oráculo no hay columnas numéricas comparables (`sin datos` en `10_equivalencia`). |
| `03_f01_bioactividades_por_tag`<br>bioactividades_por_tag() | `03_calidad_bioactividades.ipynb` | `03_bioactividades_por_tag.csv` | — | Contra el oráculo v4: `cluster_consistent` idéntico; los conteos cambian como explica el recorte. |
| `03_f02_promiscuidad`<br>promiscuidad() | `03_calidad_bioactividades.ipynb` | `03_promiscuidad_por_compuesto.csv`<br>`03_curva_filtrado.csv` | — | Contra el oráculo v4: de 30 a 7 promiscuos, explicado por el recorte. La corrida exige el mapa de compuestos para no producir la lista falsa (PROVENANCE §3.10). |

## genome_prioritization

Funciones en `genome_prioritization/figuras_genome.py`.

| figura | notebook | lee | aparece en | verificación |
|---|---|---|---|---|
| `01_f01_positivos_por_especie`<br>positivos_por_especie() | `01_barrido_parametros.ipynb` | `01_resumen_especies.csv` | — | Deriva de la base; el total de druggables lo fija `test_resultados`. |
| `01_f02_control_v5`<br>control_v5() | `01_barrido_parametros.ipynb` | `01_control_v5.csv` | — | Es el control externo (CONVENCIONES §3): `test_los_optimos_reproducen_el_control` exige los 16 óptimos iguales y ΔAUC01 < 1e-9. |
| `01_f03_auc01_por_especie`<br>auc01_por_especie() | `01_barrido_parametros.ipynb` | `01_[0-9]*.csv` | informe, presentación, página | Barrido idéntico al oráculo v4 (149 de 149 columnas); óptimos iguales al control (`test_resultados`); `test_fuga` (sin fuga en el leave-one-species-out) y `test_metricas` (pAUC y McClish). |
| `02_f01_roc_26`<br>roc_foco() | `02_optimo_y_consistencia.ipynb` | `02_ranking_26.csv`<br>`02_optimos_por_especie.csv` | informe, página | La AUC01 de la leyenda se recalcula del ranking guardado y coincide con la del óptimo en `02_optimos_por_especie`; `test_metricas`, `test_fuga`. |
| `02_f02_grilla_beta`<br>grilla_beta() | `02_optimo_y_consistencia.ipynb` | `01_[0-9]*.csv` | informe, página | El barrido que resume es idéntico al del oráculo v4 (`10_equivalencia`). |
| `02_f03_plateau`<br>plateau() | `02_optimo_y_consistencia.ipynb` | `02_plateau.csv` | informe, página | `02_plateau` idéntica a la del oráculo v4. |
| `02_f04_consistencia_especies`<br>consistencia_especies() | `02_optimo_y_consistencia.ipynb` | `02_spearman_especies.csv` | página | `02_spearman_especies` idéntica a la del oráculo v4. |
| `02_f05_enriquecimiento_topk`<br>enriquecimiento_topk() | `02_optimo_y_consistencia.ipynb` | `02_consistencia.csv` | — | `02_consistencia` idéntica a la del oráculo v4. |

## huerfanas

Funciones en `huerfanas/figuras_huerfanas.py`.

| figura | notebook | lee | aparece en | verificación |
|---|---|---|---|---|
| `01_f01_grupos_por_especie`<br>grupos_por_especie() | `01_pseudohuerfanas.ipynb` | `01_k1_drugs.csv` | — | Sin control directo de `01_k1_drugs`; el conjunto evaluado tiene el mismo tamaño que en el oráculo v4. |
| `01_f02_frank_por_grupo`<br>frank_por_grupo() | `01_pseudohuerfanas.ipynb` | `01_*_pseudohuerfanas.csv`<br>`01_k1_drugs.csv` | — | `test_fuga`: la orfanización quita todas las aristas de la droga. Contra el oráculo v4 cambia en la dirección que predice el recorte. |
| `01_f03_semilla_por_especie`<br>semilla_por_especie() | `01_pseudohuerfanas.ipynb` | `01_*_pseudohuerfanas.csv` | informe, presentación, página | Consistencia interna: con semilla no informativa el puntaje del blanco es 0 por construcción, y se observa 0 % recuperado. Contra el oráculo v4 cambia como predice el recorte. |
| `01_f04_recuperacion_por_semilla`<br>recuperacion_por_semilla() | `01_pseudohuerfanas.ipynb` | `01_resumen_frank.csv` | informe, página | La misma consistencia interna que `01_f03`, sobre `01_resumen_frank` (comparada contra el oráculo v4). |
| `02_f02_palancas`<br>palancas() | `02_cobertura_semilla.ipynb` | `02_palancas.csv` | — | Sin control. Dos de las tres palancas no pueden rescatar nada en esta base por construcción, así que el 0 es esperado (ver informe). |
| `02_f03_techo_teorico`<br>techo_teorico() | `02_cobertura_semilla.ipynb` | `02_cobertura_por_especie.csv` | — | Es la fracción informativa de `huerfanas/01` por organismo: coincide con `01_f03` (consistencia interna). |
| `03_f01_distribucion_rg`<br>distribucion_rg() | `03_desorfanizacion_global.ipynb` | `03_ranking_global.csv` | — | Contra el oráculo v4 (`03_ranking_global`): cambia como predice el recorte. |
| `03_f02_recuperacion`<br>recuperacion() | `03_desorfanizacion_global.ipynb` | `03_curva_recuperacion.csv`<br>`03_rg_estrella.csv` | informe, página | Contra el oráculo v4: r*G cambia como predice el recorte. No hay control independiente del corte. |
| `03_f03_rss`<br>rss() | `03_desorfanizacion_global.ipynb` | `03_rss_distribucion.csv` | — | Contra el oráculo v4 (`03_rss_distribucion`): cambia como predice el recorte. |
| `03_f04_directa_indirecta`<br>directa_indirecta() | `03_desorfanizacion_global.ipynb` | `03_ranking_global.csv` | — | Contra el oráculo v4 (`03_clases_inferencia`). La partición 68/32 del trabajo de 2016 no se verificó contra el artículo. |
| `03_f05_semilla_informativa`<br>semilla_informativa() | `03_desorfanizacion_global.ipynb` | `03_informativa.csv`<br>`03_informativa_topk.csv` | informe, presentación, página | Sin control independiente. Referencias internas: el azar en el panel (a); el efecto del tamaño de semilla se comprobó dentro de cada organismo, no sólo en el total. |
| `04_f01_embudo`<br>embudo() | `04_aplicacion_especie.ipynb` | `04_embudo.csv` | — | Contra el oráculo v4 (`04_sp26_embudo`): menos tratables, explicado por el recorte. |
| `04_f02_rg_sugerencias`<br>rg_sugerencias() | `04_aplicacion_especie.ipynb` | `04_sugerencias.csv`<br>`03_rg_estrella.csv` | — | Sin sugerencias bajo r*G: no se dibuja. |
| `04_f03_familias`<br>familias() | `04_aplicacion_especie.ipynb` | `04_familias.csv` | — | Sin sugerencias bajo r*G: no se dibuja. |

## verificacion

Funciones en `verificacion/figuras_verificacion.py`.

| figura | notebook | lee | aparece en | verificación |
|---|---|---|---|---|
| `10_f01_equivalencia`<br>equivalencia() | `10_equivalencia.ipynb` | `10_equivalencia.csv` | informe, presentación, página | Es la verificación: compara cada tabla contra `oraculo_v4/` con tolerancia relativa 1e-9 (`verificacion/funciones_verificacion.py`). |
