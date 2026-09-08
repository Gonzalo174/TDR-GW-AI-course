# Plan de trabajo — proyecto final GW-AI

> **Estado al 2026-09-08.** Las fases 1 y 2 están hechas: `comun/tdr.py` corre
> contra la base codificada, las 35 pruebas de `comun/tests/` pasan y
> `resultados/` está creada. Lo que sigue es la fase 3, correr los notebooks.

Este documento tiene dos partes. La **A** es el inventario de lo que hay que
tocar para que los análisis de `TDR_2026_v4` corran sobre la DB codificada de
este repositorio. La **B** es el plan para completar el proyecto final según
`final-project.html`.

Todo número que aparece acá fue medido, no recordado; el comando que lo produjo
está indicado al lado.

---

## Parte A — Qué hay que modificar

### A.0 Lo que se copió

De `TDR_2026_v4`, excluyendo `createDB/` (la base no se regenera acá):

| carpeta | contenido | LOC de funciones |
|---|---|---|
| `comun/` | `tdr.py` (rutas, carga, métricas, figuras), `nucleo.py` (el modelo, congelado), `tests/` (5 suites) | 1053 |
| `analiceDB/` | 3 notebooks: dimensiones de la red, capas y conectividad, calidad de bioactividades | 331 |
| `genome_prioritization/` | 2 notebooks: barrido de parámetros, óptimo y consistencia | 374 |
| `huerfanas/` | 4 notebooks: pseudohuérfanas, cobertura de semilla, desorfanización global, aplicación a *P. falciparum* | 562 |
| `paper/` | `informe.tex`, `refs.bib`, figuras y tablas ya generadas por v4 | — |

No se copió `bibliografia/` (PDFs de artículos de terceros) ni los diccionarios
de `acondicionarDB/mapeos/`.

**Las salidas de v4 en `paper/datos/` son el oráculo del port.** Son 40 tablas
CSV producidas por estos mismos notebooks contra la base sin codificar. Cualquier
resultado nuevo se compara contra ellas; ahí está la evidencia que pide la
consigna.

### A.1 El problema central: la DB nueva es toda enteros

    VIEJA  target_id=A0A504X1N1  sp_id=ldon  activity_tag=positive  type=Domain
    NUEVA  target_id=89920       sp_id=19    activity_tag=2         type=5

Esto no rompe con una excepción: **rompe en silencio**. Medido sobre
`DB/03_bioactivities_target_compound.csv` con pandas 2.2.3:

```
b["activity_tag"] == "positive"   ->        0 filas   (el patrón de cargar_db)
(b["activity_tag"] == 2) & cc     ->  248 457 filas   (lo correcto)
```

Un `cargar_db()` sin tocar devuelve `posdt` vacío, y todo el pipeline sigue
adelante produciendo NaN y rankings vacíos. Es el riesgo número uno del port y
la razón por la que el paso A.6 (tests que fallen antes de que fallen los
notebooks) va primero y no último.

### A.2 Inventario de literales a traducir

| # | qué | dónde | archivos |
|---|---|---|---|
| 1 | `"positive"` / `"negative"` / `"indeterminate"` / `"inconsistent"` | filtros de bioactividad | `comun/tdr.py`, `huerfanas/funciones_huerfanas.py`, `comun/tests/test_integridad.py`, `analiceDB/03_*.ipynb` |
| 2 | `type == "Domain"` (y `Homologous_superfamily`) | filtro de InterPro en `cargar_db` | `comun/tdr.py` |
| 3 | códigos de especie `'pfal'`, `'hsap'`, `'tcr'`… (41 apariciones) | `SPECIES`, `NOMBRE_CORTO`, selección de especie en notebooks | `comun/tdr.py` + 6 notebooks |
| 4 | `ann != "-1"` (faltante de OrthoMCL) | `cargar_db` | `comun/tdr.py` |
| 5 | `sp_id != "bioactive"` | `especies_con_druggables` | `comun/tdr.py` |

**Resuelto (2026-09-08).** No se publicó ningún diccionario. Las equivalencias
que el modelo necesita —y solo esas— quedaron escritas a mano en `comun/tdr.py`,
cada una con el comentario de qué selecciona:

```python
TAG_POSITIVE = 2       # se seleccionan las bioactividades positivas
IPR_DOMAIN   = 5       # se selecciona el tipo Domain
SP_BIOACTIVE = 14      # pseudo-especie, fuera de toda cuenta por especie
FALTANTE     = -1
ESCALA_PESO  = 100     # los pesos van escalados x100
```

Las especies se identifican por código. Lo que sí se publica es el **tipo de
organismo** (grupo, reino, parásito), en `SPECIES`, porque el modelo lo usa y sin
él las figuras por grupo no se leen. `SP_FOCO = 26` nombra por código el
protozoo parásito de `huerfanas/04`.

Los diccionarios grandes —`mapa_blanco` (UniProt), `mapa_compuesto`,
`mapa_cluster_*`, `mapa_interpro`, `mapa_orthomcl`— siguen privados.

**Fuga aceptada:** tres de las 16 especies son la única de su combinación
(grupo, parásito), así que para ellas el par identifica al organismo. Es la
consecuencia de publicar el grupo, y está anotada en el código.

### A.3 Verificado: la DB nueva conserva las 16 especies del modelo

```
30 especies en 00_specie_target (la base trae 30; el modelo usa 16)
17 con más de 10 druggables  =  las 16 del modelo + el pseudo-sp "bioactive" (código 14)
```

`especies_con_druggables()` excluye `bioactive` por nombre, así que con enteros
ese filtro deja de funcionar y devolvería 17. Es el punto 5 de la tabla A.2.
La buena noticia es que **las 16 especies sobreviven al recorte**: el conjunto
del modelo no cambia y las comparaciones contra `paper/datos/` son legítimas.

### A.4 Rutas: cinco constantes fuera del repo

En `comun/tdr.py`:

```python
V4         = /home/ggiordano/TDR/TDR_2026_v4      # raíz del proyecto
DB         = /data1/TDR-2025/TDR_v7/gon3/DB       # base real, no la de este repo
RAW        = /data1/TDR-2025/TDR_v7/raw_data      # .tsv de InterProScan
GON4       = V4/gon4 -> /data1/.../gon4           # salidas pesadas
CONTROL_V5 = /data1/.../gon3/resultados/genoma_completo_v5
```

Ninguna existe para quien clone el repositorio, y la consigna pide exactamente
eso ("readable by something that has never spoken to you"). Hay que:

- derivar todo de la raíz del repo (`Path(__file__).resolve().parents[1]`);
- apuntar `DB` al `DB/` de este repo;
- mover las salidas a `resultados/` dentro del repo, ignorado por git salvo las
  tablas chicas que alimentan el paper;
- `RAW` solo lo usa `nombres_ipr()`, que ya tiene cacheado
  `comun/datos_derivados/ipr_nombres.csv`. Con la DB codificada esos nombres no
  se pueden ligar a un `ann` entero sin `mapa_interpro`, que es privado: la
  función queda inutilizable y hay que decidir si se elimina o si se publica una
  versión del cache indexada por código;
- `CONTROL_V5` pesa **196 KB en 16 CSV**: entra en el repo sin problema y vale
  la pena, porque es un control independiente de los óptimos por especie.

### A.5 La capa química cambió de forma

`cargar_db(quimica=True)` hace `pd.read_csv(DB/"01_edges_clusters_fingerprint.csv")`,
que ya no existe: ahora son dos `.gz` que hay que concatenar. Y el peso pasó a
entero ×100 (`0.82 -> 82`), así que todo umbral sobre `weight` hay que
reescalarlo. La usan 4 notebooks (`analiceDB/02`, `huerfanas/01`, `02`, `04`).

### A.6 `cluster_consistent` dejó de ser booleano

Era `True`/`False`, ahora es `1`/`0` int64. En pandas 2.2.3 el `pos &= cc` de
`cargar_db` **no** falla con enteros, lo cual es peor que si fallara: conviene
castear a bool explícitamente en la carga.

### A.8 Los dos vocabularios de anotación se pisan

Lo encontró `comun/tests/test_integridad.py` al correrlo por primera vez contra
la base nueva: 6 enlaces proteína-categoría duplicados donde no debía haber
ninguno. La causa resultó más grave que el síntoma.

InterPro y OrthoMCL se codificaron **por separado, cada uno arrancando en 0**.
Medido: los 14 083 códigos de dominio InterPro caen todos dentro del rango de
OrthoMCL (0–76 156). El entero 100 es a la vez un dominio y un grupo de
ortología que no tienen nada que ver, y `cargar_db` los concatenaba en la misma
columna `ann`.

El daño no es sólo el duplicado. `nucleo.get_annot_druggability_pv` corre el
test de Fisher **por vocabulario**, y elige las filas de cada uno con
`ann.str.startswith("IP")` y `("OG")`. Sobre enteros eso no selecciona nada, y
sobre códigos fusionados el conteo de proteínas por categoría queda mal de raíz.

**Arreglado** reponiendo el prefijo en la carga (`"IP" + código`, `"OG" + código`).
Restaura la semántica original exactamente y deja `nucleo.py` sin tocar. Hay un
test nuevo, `test_los_dos_vocabularios_no_se_pisan`, que lo vigila.

### A.9 Un cambio a `nucleo.py`, el módulo congelado

`get_druggable_targets` arma la clave de su diccionario de salida con
`"_".join(sp_out)`, que sobre códigos enteros lanza `TypeError`. Se cambió por
`"_".join(map(str, sp_out))`: arma un nombre, no interviene en ningún cálculo.
Es la única línea modificada del modelo y está comentada como tal.

El envoltorio `tdr.semilla_global` absorbe el resto: `nucleo` envolvía `sp_out`
en lista sólo si era `str`, cosa que dejó de cumplirse.

### A.7 Lo que cambia de resultado, y por qué está bien

El recorte de la base (componentes conexas sin ningún cluster con bioactividad
positiva: 2 396 103 → 831 175 compuestos) **cambia legítimamente** los
descriptivos de `analiceDB/`, que justamente miden esas componentes. No así
`genome_prioritization/` ni `huerfanas/`, que trabajan sobre la capa de
anotaciones y los druggables: ahí los números **deberían reproducirse**.

Esa asimetría es el checkeo más fuerte que tenemos y hay que explotarla:
lo que debe cambiar, cambia de forma explicable; lo que no debe cambiar, no cambia.

---

## Parte B — Plan para el proyecto final

La consigna pide cuatro cosas: **provenance** de cada resultado, **checks**
detrás de cada figura, una **página HTML** y un **PDF** para presentar, y un
repositorio que **un agente pueda reproducir sin haber hablado con nadie**.

### Fase 1 — Reproducibilidad de la infraestructura ✅

1. `comun/rutas.py`: raíz derivada del archivo, cero paths absolutos.
2. `DB/mapeos/` con los cinco diccionarios chicos + `comun/codigos.py` que
   traduce en un solo lugar (A.2).
3. `cargar_db()` adaptado: tags, tipos, `-1`, `bioactive`, `cluster_consistent`
   a bool, y `leer_aristas()` que concatena los `.gz` y reescala el peso (A.5).
4. `requirements.txt` con versiones exactas (hoy: pandas 2.2.3, python 3.x, y
   la restricción de glibc 2.27 que impide pyarrow >15, ya documentada).
5. Copiar `CONTROL_V5` (196 KB) a `control/`.

**Hecho.** `RAIZ = Path(__file__).resolve().parents[1]` y todo cuelga de ahí;
`DB` apunta al `DB/` del repo; `out()` escribe en `resultados/<carpeta>_out/`,
ya creada y versionada; `leer_aristas()` concatena los `.gz` y desescala el peso;
`RAW` quedó como ruta opcional por `TDR_RAW` (los crudos no se publican);
`requirements.txt` fija las versiones; `control/genoma_completo_v5/` copiado con
los archivos renombrados a código (venían como `pfal.csv`).

Verificado de punta a punta: `cargar_db()` da **248 457** positivos y **101 885**
negativos, `especies_con_druggables()` devuelve exactamente **16**, y
`leer_aristas()` reconstruye **36 152 622** aristas con `weight` en [0.81, 0.99].

### Fase 2 — Tests que atrapen la falla silenciosa ✅

Adaptar las 5 suites de `comun/tests/` y agregar las que faltan para el port:

- `posdt` no vacío y con el conteo exacto **248 457** (medido en A.1);
- 16 especies con druggables, sin `bioactive`;
- la tabla de aristas concatenada tiene **36 152 622** filas (`DB/meta.json`);
- toda columna de `DB/` es entera (invariante que declara `acondicionarDB`);
- las métricas (`pauc_normalizada`, `mcclish`) contra valores de referencia.

**Hecho.** `python3 -m unittest discover -s comun/tests -p "test_*.py"`: **35
tests, todos en verde** (2 salteados: la capa química, que se activa con
`TDR_TESTS_QUIMICA=1`). Cuatro son nuevos y cubren justamente lo que el port
puede romper en silencio: `test_el_filtro_de_positivos_no_quedo_vacio`,
`test_la_base_es_entera`, `test_los_dos_vocabularios_no_se_pisan` y
`test_la_base_esta_y_esta_completa`.

La primera corrida encontró 6 fallas reales, entre ellas la de A.8.

### Fase 3 — Correr los 9 notebooks y comparar contra el oráculo

En orden: `analiceDB` (3), `genome_prioritization` (2), `huerfanas` (4).
Cada uno escribe su `NN_meta.json` (ya está en el diseño de v4).

Después, un notebook nuevo — `10_equivalencia.ipynb` — que compare tabla por
tabla contra `paper/datos/` y clasifique cada diferencia en:

| clase | qué significa |
|---|---|
| idéntico | el port no alteró nada |
| explicado por el recorte | `analiceDB/`, con el número de compuestos perdidos que lo justifica |
| **discrepancia** | hay que investigarlo antes de seguir |

**Entregable:** la tabla de equivalencia. Es la pieza que sostiene el
"defenderlo ante alguien que no te vio producirlo".

### Fase 4 — Provenance explícito

Por cada figura y cada número del informe: qué notebook y qué función lo
produjo, qué entró, qué es código propio y qué es librería, y qué alternativa
defendible se descartó. Parte ya existe (`meta.json` por corrida, README §1.4);
falta el índice que lo recorra completo y el registro de las decisiones
—recorte, codificación, empaquetado— que ya están escritas en
`acondicionarDB/README.md` y hay que traer al repo.

### Fase 5 — Los dos entregables para personas

- **PDF:** `paper/informe.tex` ya compila; hay que actualizar los números que
  cambien por el recorte y agregar la sección de reproducibilidad.
- **HTML:** una página que cuente el problema, el método, las figuras y el
  camino de verificación, con los enlaces al repo. Es la que se va a linkear
  desde la página general del curso.

### Fase 6 — Prueba de fuego

Clonar el repositorio en un directorio limpio y reproducirlo siguiendo solo el
`README.md`, sin usar nada de la máquina original. Lo que falle ahí es lo que le
va a fallar a quien lo corrija.

### Orden y dependencias

```
Fase 1 ──> Fase 2 ──> Fase 3 ──> Fase 4 ──> Fase 5 ──> Fase 6
              │            │
              └── rojo ────┘   (un test que falla vuelve a Fase 1)
```

### Decisiones abiertas

1. **`nombres_ipr()`**: sin `mapa_interpro` no hay descripción legible de los
   dominios. ¿Se elimina la función, o se publica el cache reindexado por código?
   Afecta a `huerfanas/04`, que hoy nombra los dominios de las familias que
   propone. Por ahora la función quedó con la advertencia en el docstring.
4. **`paper/`**: el informe de v4 nombra las especies en el texto, en las figuras
   y en los nombres de archivo (`01_pfal.csv`, `01_calb_pseudohuerfanas.csv`).
   Si el criterio de no identificar especies vale para todo el repositorio, hay
   que renombrarlos y reescribir el `.tex`; si vale sólo para el código y los
   datos, `paper/` queda como está.
2. **Alcance**: ¿el proyecto final son los tres análisis completos, o uno solo
   hecho a fondo? La consigna valora "many components" pero también pide
   evidencia suficiente para defender cada resultado.
3. **`paper/informe.tex`** está escrito contra los números de v4. ¿Se actualiza,
   o el informe del curso es un documento nuevo y más corto?
