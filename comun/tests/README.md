# Pruebas de `comun/tests/`

Escritas con `unittest` de la biblioteca estándar. `nucleo.py` está congelado
(CONVENCIONES.md §1): se cambia sólo con una prueba de esta carpeta que lo
justifique. Desde la raíz del repositorio:

```bash
python -m unittest discover -s comun/tests -p "test_*.py"             # todas, ~25 s
python -m unittest discover -s comun/tests -p "test_metricas.py" -v   # un archivo
TDR_TESTS_QUIMICA=1 python -m unittest discover -s comun/tests -p "test_integridad.py"
                                                                      # + capa química (36 M aristas)
```

o `make pruebas`.

| archivo | qué verifica | tarda |
|---|---|---|
| `test_metricas.py` | McClish: azar → 0.5, perfecto → 1.0, monotonía (el ranking de parámetros no cambia), que la pAUC sólo dependa de FPR ≤ 0.1, y deja documentado por qué la fórmula de la versión v3 era incorrecta | instantáneo |
| `test_proyeccion.py` | Que `nucleo.nds()` sea equivalente a la proyección `M S Mᵀ`: con λ = 0, `score_i · k_i` coincide exactamente con `(M S Mᵀ w)_i`. También que β = 1 penalice a las categorías grandes | instantáneo |
| `test_rutas.py` | Que todo cuelgue de la raíz del repositorio, que la base esté completa, que `out()` rechace un nombre de carpeta con error, que ninguna carpeta arranque con dígito y que todo notebook lleve su `NN_` | instantáneo |
| `test_integridad.py` | Tablas de entrada: la base es entera, el filtro de positivos no queda vacío, toda proteína anotada tiene especie, los prefijos `IP`/`OG` no se pisan, sin duplicados, positivos y negativos disjuntos, pesos en rango | ~10 s (+ química, minutos) |
| `test_vecinos.py` | Que la semilla de una droga (`use_neighbors=True`) busque vecinos en `ddt`/`dds` con el `cluster_id` de la droga y no con su `drug_id`: con datos sintéticos, trae el vecino real y no el del cluster que tiene el mismo número que la droga. Falla con el código anterior al 2026-09-24 (PROVENANCE.md §3.16) | instantáneo |
| `test_fuga.py` | Que `sp_out` deje la semilla sin ningún blanco de la especie evaluada (con control negativo) y que la orfanización remueva todas las aristas de la droga | ~5 s |
| `test_entorno.py` | Que todo lo que el repositorio importa, notebooks incluidos, esté instalado y fijado en `requirements.txt`, y que toda celda parsee | ~2 s |
| `test_resultados.py` | Sobre las tablas de `resultados/`: los 16 óptimos reproducen el control, todas las especies superan el azar, y el recorte no movió la capa de anotaciones ni la evidencia positiva | instantáneo |
| `test_figuras.py` | Que toda figura se dibuje desde sus tablas con la carga de la base anulada, que no lea lo que no declaró, que declare qué la verifica, que `FIGURAS.md` esté al día, que ningún notebook dibuje por su cuenta y que los entregables sólo usen figuras registradas | ~10 s |

Las que leen la base se saltean solas si `DB/` no está, y las de resultados, si
faltan las tablas (están versionadas, así que en un clon corren).

**Orden de trabajo** (CONVENCIONES.md §9): estas pruebas van en verde antes de
correr cualquier notebook. De las cuatro fallas silenciosas del port, dos las
encontraron estas pruebas y dos la inspección de las salidas (PROVENANCE.md §3).
