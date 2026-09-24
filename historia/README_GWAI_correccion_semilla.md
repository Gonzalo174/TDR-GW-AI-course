> **Documento histórico — aplicado el 2026-09-24.** Lo escribió una sesión de
> Claude que trabajaba en `TDR_2026_v4` (la versión sin codificar del proyecto),
> para que otra sesión portara a este repositorio la corrección que había hecho
> allá. Se conserva tal cual, debajo de esta nota, como registro de esa
> interacción entre las dos versiones. Sus rutas (`/home/ggiordano/TDR/...`,
> `gon4/`) no son de este repositorio.
>
> **Qué se hizo acá, siguiendo el §5:**
>
> 1. La corrección de `comun/nucleo.py` (§2), con `comun/tests/test_vecinos.py`
>    copiado de v4. La prueba **falla con el código anterior** (2 fallas y 1
>    error, igual que en v4) y pasa con el nuevo. Antes se confirmó sobre
>    `DB/` que `ddt` y `dds` usan ids de cluster y que las numeraciones se
>    solapan. Registro: `PROVENANCE.md` §3.16 y `CONVENCIONES.md` §1.
> 2. **Destino de los resultados:** se reemplazó `resultados/huerfanas_out/` en
>    lugar de crear una carpeta nueva, porque el informe, las figuras y
>    `make verificar` leen de ahí. La corrida anterior no se perdió: está en el
>    commit `e3afdb0`.
> 3. Se re-corrieron `huerfanas/01`–`04` y `verificacion/10`. `genome_prioritization/`
>    no se tocó. `analiceDB/03` se re-corrió por la lista de promiscuos (punto 5).
> 4. **La predicción del §3 se cumple.** Semilla nula 67.8 % (v4 corregida:
>    67.2 %), 1 365 semillas informativas (1 383), recuperadas 17.5 % (17.7 %),
>    r*G 76 (80), recuperadas directa/indirecta 1 113 / 246 (1 124 / 250). Cruzando ids con el mapa privado:
>    en los 5 218 compuestos que comparten las dos muestras la clase de semilla
>    coincide en el 100 %, y frank en todos menos uno. Las dos sugerencias para
>    *P. falciparum* son los mismos compuestos, el mismo blanco, el
>    mismo puntaje y el mismo rG que en v4 corregida. Así se cerró el pendiente de
>    `PLAN.md` sobre la semilla nula (83 contra 69 %): la causa era este error, no
>    el filtro de promiscuidad. Detalle en `PROVENANCE.md` §5.
> 5. **§4.1, promiscuos duplicados:** confirmado (7 filas, 5 compuestos). Se
>    corrigió `derivar.py`, que ahora deduplica, y se re-corrió con los datos
>    privados. Las 7 445 relaciones filtradas no cambian, y el filtro de
>    `huerfanas/` ya usaba la lista como conjunto.
> 6. **§4.4, sensibilidad de r*G:** `huerfanas/03` escribe ahora
>    `03_rg_sensibilidad.csv`. Con k_σ entre 2 y 4, r*G va de 68 a 88. El informe
>    lo presenta como un corte operativo.
> 7. §4.3 (`desde = 50`) y §4.5 (afirmaciones sin respaldo): el texto de este
>    repositorio no las tenía. La partición 68/32 atribuida al paper ahora se
>    rotula como no verificada.
> 8. Apareció una falla nueva: `huerfanas/04` no había producido nunca una
>    sugerencia, así que el cruce con las anotaciones nunca se había ejecutado,
>    y ahora fallaba por tipo (entero contra texto). Registro: `PROVENANCE.md` §3.17.
> 9. **Ids censurados en esta copia.** El original citaba los ids sin codificar
>    de los compuestos sugeridos y los accessions de sus blancos; junto con
>    `04_sugerencias.csv` habrían revelado su código (PROVENANCE.md §3.13). Acá se
>    reemplazaron por conteos.
> 10. Se regeneraron figuras, `numeros.tex`, el informe, la presentación y la página, y se
>    revisó el texto. `make verificar` pasa: 63 pruebas (2 salteadas), 36 de 36
>    figuras, 356 columnas contra el oráculo sin discrepancias.

---

# README_GWAI — correcciones pendientes de portar desde TDR_2026_v4

> **Para la sesión de Claude que lea esto.** Este documento lo escribió otra sesión
> (2026-09-24), trabajando en `/home/ggiordano/TDR/TDR_2026_v4` (la versión sin
> codificar del mismo proyecto). Ahí se encontró y corrigió un error en
> `comun/nucleo.py` que **este repositorio también tiene** (verificado en el código,
> no supuesto). Este archivo explica el error, cómo corregirlo acá, qué se va a mover
> y cómo verificarlo. **Nada de lo que describe está aplicado todavía en este
> repositorio.** El archivo no está commiteado.
>
> **Antes de tocar nada, confirmá con el usuario:** si quiere la corrección acá, dónde
> guardar los resultados nuevos (en v4 pidió **no sobrescribir** y crear una carpeta
> nueva, ver §6) y si commitear. `resultados/` está versionado y `make verificar`
> compara el informe contra esas tablas.

---

## 1. El error: la semilla busca vecinos químicos con el id equivocado

**Dónde:** `comun/nucleo.py::get_druggable_targets`, rama `use_neighbors=True`
(sólo se usa cuando se pide la semilla de compuestos concretos, `doi=[...]`).
En este repositorio: líneas **161–192** (hoy, commit `e3afdb0`).

```python
doi_and_nei = np.concatenate([nei1, [idoi]]).astype(float)          # ids de COMPUESTO
neitc_1 = ddt[ddt["clusID1"].isin(doi_and_nei)][["clusID2", "weight"]].copy()  # ids de CLUSTER
neitc_2 = ddt[ddt["clusID2"].isin(doi_and_nei)][["clusID1", "weight"]].copy()
...
neisc = dds[dds.index.isin(doi_and_nei)][["to", "weight"]].copy()   # dds indexado por "from" (cluster)
```

- `ddt` (`01_edges_clusters_fingerprint*`) y `dds` (`02_edges_clusters_subestructure`)
  unen **clusters identitarios**: `clusID1/clusID2` y `from/to` son `cluster_id`
  de `01_clusters_fingerprint` / `02_clusters_subestructure`. Así los arma
  `createDB/01` en v4 (`edges_save = edges_agg[["A","B",...]].rename(columns={"A":"clusID1","B":"clusID2"})`,
  donde A y B son ids de cluster identitario).
- El código los consulta con `doi_and_nei`, que son **ids de compuesto** (la droga y
  sus *cluster-mates*). Viene de `nds_fun.R` (`setkey(ddt, clusID1); ddt[.(as.numeric(c(nei1, doi)))]`),
  escrito cuando esas tablas eran droga–droga. Se tradujo tal cual a Python.
- **Falla en silencio** porque las numeraciones se superponen. En este repositorio
  (base codificada): compuestos `0..831 174`, clusters de fingerprint `0..769 636`,
  clusters de subestructura `0..788 536`. Buscar el compuesto 12345 como si fuera el
  cluster 12345 devuelve los vecinos **de otro cluster**, que no tiene relación con
  la droga.
- **No afectado:** los *cluster-mates* (tipo `identity`: `ssc`, `stc`, `stc2`,
  `ssc2`, que usan `tclus`/`sclus` correctamente) y todo lo que no pasa por esta rama.

### Qué se ve afectado y qué no

| | ¿afectado? | por qué |
|---|---|---|
| `genome_prioritization/` (barrido, óptimos, β, quinasas) | **no** | la semilla viene de `tdr.semilla_global` → `get_druggable_targets(posdt, st, sp_out=...)` **sin `doi`**: entra en la rama "todas las drogas" y nunca toca `ddt`/`dds`; además esos notebooks cargan la base sin la capa química |
| `huerfanas/01` (`fh.priorizar_droga` → `fh.semilla_de_droga`) | **sí** | clases de semilla, frank, rG, rSS, directa/indirecta |
| `huerfanas/02` (`fh.cobertura_por_palanca` → `semilla_de_droga`) | **sí** | |
| `huerfanas/03` (r*G, curvas, `03_informativa*`) | **sí** | lee lo de 01 |
| `huerfanas/04` sugerencias (`04_sugerencias.csv`, `04_familias.csv`) | **sí** | el **embudo** (`fh.embudo_especie`) **no**: usa `ddt` con ids de cluster, correctamente |
| `fh.construir_pseudohuerfanas` (grupos 0–3) | no | usa `ddt` con ids de cluster |
| `verificacion/10_equivalencia` | cambia | el oráculo v4 tiene el mismo error; las columnas de `huerfanas_out` van a dejar de coincidir y esa diferencia es la corrección |
| `informe/`, `index.html`, presentación | **sí** | todo lo de desorfanización (`numeros.py` lee de `resultados/`) |

---

## 2. La corrección (la que se aplicó en v4)

Un único salto en la capa química, **como en el original**. El usuario pidió
explícitamente **no** agregar un segundo salto. Sólo cambia el id con el que se consulta:

```diff
                 doi_and_nei = np.concatenate([nei1, [idoi]]).astype(float)
+                # ddt y dds unen CLUSTERS identitarios (createDB/01-02): se consultan
+                # con los cluster_id de estas drogas, no con sus drug_id. Hasta el
+                # 2026-09-24 se consultaban con drug_id (heredado de nds_fun.R, cuando
+                # las aristas eran droga-droga) y devolvian los vecinos de otro
+                # cluster. Test: comun/tests/test_vecinos.py
+                t_cl = tclus.loc[tclus["drug"].isin(doi_and_nei), "cluster_id"].unique()
+                s_cl = sclus.loc[sclus["drug"].isin(doi_and_nei), "cluster_id"].unique()

-                neitc_1 = ddt[ddt["clusID1"].isin(doi_and_nei)][["clusID2", "weight"]].copy()
+                neitc_1 = ddt[ddt["clusID1"].isin(t_cl)][["clusID2", "weight"]].copy()
 ...
-                neitc_2 = ddt[ddt["clusID2"].isin(doi_and_nei)][["clusID1", "weight"]].copy()
+                neitc_2 = ddt[ddt["clusID2"].isin(t_cl)][["clusID1", "weight"]].copy()
 ...
-                neisc = dds[dds.index.isin(doi_and_nei)][["to", "weight"]].copy()
+                neisc = dds[dds.index.isin(s_cl)][["to", "weight"]].copy()
```

Referencia: `/home/ggiordano/TDR/TDR_2026_v4/comun/nucleo.py` (ya corregido; el
docstring del módulo lista el cambio). Las reglas del núcleo (máximo por
cluster y blanco sumado sobre clusters, `len(cdwt) > 1`, `len(a) > 1`, filtro de
evidencia contradictoria) **no se tocaron**.

### Qué exigen las convenciones de este repositorio

- `CONVENCIONES.md §1`: `nucleo.py` está congelado y sólo se cambia **con una prueba
  que lo justifique**. Hoy la única excepción es `PROVENANCE.md §3.4`; habría que
  sumar esta como `§3.x` y actualizar la frase "se cambió una línea" de §1.
- La prueba ya existe en v4 y usa datos sintéticos (no depende de `DB/`):
  `/home/ggiordano/TDR/TDR_2026_v4/comun/tests/test_vecinos.py`. Se puede copiar tal
  cual a `comun/tests/`. Arma un compuesto 7 en el cluster 100 con vecino real en el
  101, y un **cluster señuelo numerado 7** con otro vecino. En v4 se comprobó que
  **pasa con la corrección y falla con el código anterior** (2 fallas y 1 error).
  Conviene repetir esa comprobación acá. Ojo con `test_entorno`/`test_rutas` si
  validan la lista de archivos de tests.

---

## 3. Qué resultó en v4 (como referencia, no como expectativa exacta)

Misma muestra de 7782 pseudohuérfanas (hasta 1000 por especie, `RSEED = 123457`),
mismos óptimos, mismo filtro de promiscuidad; sólo cambió la semilla. Tablas:
`/home/ggiordano/TDR/TDR_2026_v4/gon4/huerfanas_out` (antes) y `.../huerfanas_out2` (después).

| | v4 antes | v4 corregido | GW-AI hoy (`resultados/huerfanas_out`) |
|---|---:|---:|---:|
| semilla nula | 83,1 % | 67,2 % | ~69 % (`PROVENANCE.md`) |
| semillas informativas | 207 | 1383 | ver `01_resumen_frank.csv` |
| recuperadas (frank < 0,1) | 2,6 % | 17,7 % | ver `01_resumen_frank.csv` |
| r*G | 54 | 80 | 74 (`informe/numeros.tex`) |
| directa / indirecta (recuperadas) | 124 / 81 | 1124 / 250 | |
| sugerencias pfal (código 26) | 1 compuesto → 2 blancos | 2 compuestos → 1 blanco | (ids codificados) |

Controles que se hicieron en v4 y conviene repetir acá:
- con el código viejo, la corrida reproducía exactamente la tabla previa
  (7782/7782 en clase de semilla y frank);
- una reimplementación independiente de la semilla (`reunion_11sep2026/funciones_reunion.py::Vecindario`)
  dio la **misma** recuperación que el núcleo corregido en *A. thaliana* (22,3 %);
- las semillas nulas "sin vecinos", "vecinos sin bioactividad" y "vecinos sin
  positivos" siguen 100 % nulas después de la corrección. Sólo cambian las que
  tenían vecinos con positivos.

### Hipótesis sobre el pendiente de `PLAN.md` (sin verificar)

`PLAN.md` deja abierto *"por qué la semilla nula bajó del 83 al 69 % respecto de v4"*.
Explicación probable: con el error, la semilla depende de **qué cluster tiene el
número del compuesto**. La codificación renumera compuestos y clusters, así que los
clusters "equivocados" son otros. El recorte deja sólo componentes con alguna
bioactividad positiva, así que un cluster al azar tiene más probabilidad de aportar
blancos. **Predicción verificable:** el recorte sólo elimina componentes sin ninguna
bioactividad positiva, que no pueden aportar blancos a una semilla. Entonces, **con la
corrección**, las semillas de GW-AI deberían ser equivalentes a las de v4 para los
mismos compuestos. Si la muestra de pseudohuérfanas es la misma, que puede no serlo
porque el orden de `k1` depende de los ids, los porcentajes deberían coincidir con la
columna "v4 corregido". Si no coinciden, hay que investigar antes de reportar.

---

## 4. Otros hallazgos de v4 que aplican acá

1. **Lista de promiscuos con duplicados.** `datos_externos/promiscuidad/compuestos_promiscuos.csv`
   tiene **7 filas pero 5 compuestos distintos** (451642 y 200190 aparecen dos veces).
   En v4 pasaba lo mismo (30 filas, 25 compuestos). Revisar dónde se cita "7" en
   `informe/`, `PROVENANCE.md` y `PLAN.md`. La cuenta de aristas afectadas en v4 ya
   estaba bien deduplicada; conviene verificarlo acá también.
2. **El filtro de promiscuidad no actúa sobre la capa del modelo.** En v4 saca los
   compuestos de `sclus` pero ninguna de las aristas de `dds`. Este repositorio ya lo
   vio (`03_impacto_en_la_base`). El control "sin filtro" de v4 se hizo con la semilla
   errónea; con la corregida no se repitió.
3. **r*G se estima sobre `ℓ ≥ desde = 50`** (`fh.r_g_estrella`), no sobre `ℓ ≥ 500`.
   El informe de v4 decía 500 por error. Revisar si el texto de acá describe la cola.
4. **Sensibilidad de r*G** (v4, semilla corregida): bootstrap IC95 66–93; sólo
   variando k_σ de 2 a 4 va de 69 a 96; con toda la grilla, de 20 a 195. Conviene
   presentarlo como elección operativa, no como parámetro identificado.
5. En el informe de v4 había afirmaciones sin respaldo en los datos que se
   corrigieron. Conviene buscarlas acá:
   - "especies filogenéticamente afines" para *C. elegans*–*T. cruzi* (ρ = 0,99);
   - "consistente con el trabajo original" sin cifras;
   - "un vector global costaría poco" sin haberlo medido;
   - el nombre de IPR017452, que es "GPCR, rhodopsin-like, 7TM" y no "sitios de unión
     a ATP". Acá los dominios no se nombran, así que quizá no aplica.

---

## 5. Pasos sugeridos (una vez que el usuario confirme)

1. Aplicar el diff del §2 a `comun/nucleo.py` y copiar `test_vecinos.py`. Comprobar que
   el test **falla** con el código viejo y **pasa** con el nuevo, y correr `make verificar`
   (las pruebas que comparan contra `resultados/huerfanas_out` van a fallar: es esperado).
2. Registrar el cambio en `PROVENANCE.md` (§3.x nuevo) y en `CONVENCIONES.md §1`.
3. Re-correr `huerfanas/01`, `02`, `03`, `04` en el orden del `Makefile`, hacia la
   carpeta que el usuario elija (§6). `genome_prioritization/` y `analiceDB/` no hace
   falta re-correrlos.
4. Re-correr `verificacion/10_equivalencia` y documentar que las diferencias en
   desorfanización son la corrección y no un problema del port.
5. Regenerar figuras, `informe/numeros.tex`, informe, presentación y página, y revisar
   el **texto**: cifras escritas en prosa, interpretaciones como "el cuello de botella
   es la cobertura química" (sigue siendo cierto, pero con 67 % de nulas y no 83 %) y
   la sugerencia de la especie foco.

---

## 6. Preferencias del usuario (expresadas en la sesión de v4)

- **No sobrescribir resultados:** en v4 la corrida corregida fue a una carpeta nueva
  (`gon4/huerfanas_out2`) y la anterior quedó intacta. Acá `resultados/` está en git:
  preguntar si prefiere `resultados/huerfanas_out2/` o un commit que reemplace.
- **Un único salto** en la capa química (se evaluó un segundo salto y se descartó:
  en *A. thaliana* pasó de 22,3 % a 23,1 %).
- **Las cifras de los entregables salen sólo de las corridas.** Lo que venga de papers
  se rotula "reportado por Berenstein et al. 2016" con su ubicación en el artículo. Lo
  que no esté verificado se marca como hipótesis.
- No commitear ni pushear sin que lo pida.

## 7. Dónde está todo en v4

| qué | ruta |
|---|---|
| núcleo corregido | `/home/ggiordano/TDR/TDR_2026_v4/comun/nucleo.py` |
| test | `/home/ggiordano/TDR/TDR_2026_v4/comun/tests/test_vecinos.py` |
| corrida anterior / corregida | `/home/ggiordano/TDR/TDR_2026_v4/gon4/huerfanas_out` / `huerfanas_out2` (con `README.md` y los notebooks ejecutados en `notebooks/`) |
| controles de la reunión del 11/9 (incluye la semilla reimplementada, motivos de nulidad, sensibilidad de r*G) | `/home/ggiordano/TDR/TDR_2026_v4/reunion_11sep2026/` y `gon4/reunion_11sep2026_out/` |
| informe actualizado (sección 2.5 describe el error) | `/home/ggiordano/TDR/TDR_2026_v4/paper/informe.tex` / `.pdf` |
