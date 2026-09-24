# Insumo externo: la lista de compuestos promiscuos

Es el único cálculo del proyecto que no sale de `DB/`. El filtro de promiscuidad
química (CONVENCIONES.md §10) retira los compuestos chicos y promiscuos —peso
molecular < 150 Da y contenidos en más de 100 superestructuras, el criterio del
trabajo de 2016— antes de armar las semillas de `huerfanas/`. Contar
superestructuras exige las relaciones de subestructura **crudas**, que `DB/` no
conserva: su capa de subestructuras está clusterizada.

Por eso la lista se calcula aparte, con `derivar.py`, y se versiona acá como un
insumo más, igual que `DB/`. `analiceDB/03` la lee, la transcribe a `resultados/`
y mide desde la base cuánto saca de la capa que usa el modelo.

| archivo | qué es |
|---|---|
| `compuestos_promiscuos.csv` | los 5 compuestos que cumplen el criterio: código, superestructuras que los contienen, peso molecular |
| `impacto_filtro.csv` | cuántos compuestos y cuántas relaciones crudas saca el filtro |
| `curva_filtrado.csv` | cuántos compuestos y relaciones sacaría cada combinación de umbrales: la sensibilidad del criterio |
| `derivar.py` | el cálculo, desde los datos crudos |

## Cómo se produjo

Con `derivar.py` el 2026-09-10, sobre tres archivos privados:

| archivo | qué tiene | tamaño |
|---|---|---|
| `raw_data/subestructures_chembl35_biolip.txt` | 1 065 346 relaciones superestructura → subestructura, con ids de ChEMBL y SMILES | 248 MB |
| `raw_data/compounds/compound_data.csv` | id interno ↔ ChEMBL, SMILES y peso molecular | 293 MB |
| `acondicionarDB/mapeos/mapa_compuesto.csv` | el diccionario de id interno a código de `DB/` | 12 MB |

```bash
TDR_RAW=<carpeta raw_data> TDR_MAPEOS=<carpeta de mapeos> python derivar.py   # ~45 s
```

Los dos primeros identifican los compuestos y el tercero es lo que la
codificación oculta: ninguno se publica.

## Cómo se verificó

- `derivar.py` reproduce las tres tablas **byte a byte** respecto de las que
  había escrito `analiceDB/03` cuando todavía leía los datos crudos.
- La corrida exige el mapa: sin él, los ids crudos no cruzan contra `DB/` y la
  lista sale con compuestos falsos en vez de los 5 verdaderos, sin que nada
  falle (PROVENANCE.md §3.10).
- Contra el oráculo v4, la lista pasa de 25 a 5 compuestos, un cambio que
  explica el recorte de la base (`verificacion/10`).
- **Corrección del 2026-09-24.** Hasta esa fecha la tabla tenía 7 filas para 5
  compuestos: `compound_data.csv` repite algunos compuestos y el cruce con el
  peso molecular los duplicaba (en v4 pasaba lo mismo: 30 filas, 25
  compuestos). `derivar.py` ahora deduplica antes del cruce. El conteo de
  relaciones que saca el filtro (7 445, 1.6 %) ya estaba bien, y el filtro de
  `huerfanas/` usa la lista como conjunto, así que ninguna semilla cambió por
  esto; sí cambian el conteo de compuestos y la curva de sensibilidad.

## Qué no se publica

La tabla por compuesto (peso molecular y superestructuras de los ~84 000
compuestos): el peso molecular con toda su precisión, al lado del código,
permitiría reidentificar los compuestos y deshacer la codificación
(PROVENANCE.md §3.13). Se publican sólo los 5 filtrados y los agregados.
