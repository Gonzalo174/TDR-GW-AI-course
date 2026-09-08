# Tests de `comun/`

`nucleo.py` está congelado (README §1.1): se cambia sólo con un test de esta
carpeta que lo justifique. Están escritos con `unittest` de la biblioteca
estándar, no hace falta instalar nada.

```bash
cd /home/ggiordano/TDR/TDR_2026_v4/comun
python3 -m unittest discover -s tests -v            # todo
python3 -m unittest tests.test_metricas -v          # sólo métricas (instantáneo)
python3 -m unittest tests.test_proyeccion -v        # sólo el modelo (instantáneo)
TDR_TESTS_QUIMICA=1 python3 -m unittest tests.test_integridad -v   # + capa química (941 MB)
```

| Archivo | Qué verifica | Tarda |
|---|---|---|
| `test_metricas.py` | McClish: azar → 0.5, perfecto → 1.0, monotonía (el ranking de parámetros no cambia), que la pAUC sólo dependa de FPR ≤ 0.1, y deja documentado por qué la fórmula de `full_genome_v4.ipynb` era incorrecta. | instantáneo |
| `test_proyeccion.py` | Que `nucleo.nds()` sea equivalente a la proyección `M^bip S (M^bip)ᵀ` del roadmap: con λ = 0, `score_i · k_i` coincide exactamente con `(M S Mᵀ w)_i`. También que β = 1 penalice a las categorías grandes. | instantáneo |
| `test_rutas.py` | Que `gon4` resuelva a `/data1`, que `out()` rechace un nombre de carpeta con typo, que ninguna carpeta arranque con dígito (importabilidad) y que todo notebook lleve su `NN_`. | instantáneo |
| `test_integridad.py` | Tablas de entrada: toda proteína anotada tiene especie, prefijos `IPR`/`OG` coherentes con la columna `db` (el Fisher de `nucleo` separa por prefijo y deja NaN silenciosos si aparece otro), sin duplicados, positivos y negativos disjuntos, pesos de similitud en rango. | minutos |
| `test_fuga.py` | Que `sp_out` deje la semilla sin ningún blanco de la especie evaluada (con control negativo) y que la orfanización remueva todas las aristas de la droga. | minutos |

Los tests que leen la base se saltean solos si `DB/` no está montado. Los mismos
chequeos de integridad, hechos sobre las 8 tablas y con salida en tabla, están
en `createDB/05_qc_base.ipynb`.

**Orden de trabajo** (README §6): estos tests van en verde antes de tocar
cualquier carpeta de análisis.
