# Procedencia de `tdr.ANN_QUINASA`

Decidir si un dominio es de quinasa exige su **nombre**, y la columna `ann` de
`DB/` es un entero codificado: ligar uno con otro necesita `mapa_interpro.csv`,
que es privado. Es la limitación que `tdr.nombres_ipr()` documenta desde el
principio (PROVENANCE.md §3.6).

A diferencia de la lista de promiscuos, acá **no hay insumo versionado**: la
lista es chica y estable, así que vive directamente en `comun/tdr.py` como
`ANN_QUINASA`, con el accession de InterPro al lado de cada código. Esta carpeta
es sólo la procedencia: `derivar.py` regenera esa constante y la imprime lista
para pegar, para que la marca se pueda auditar y rehacer.

Ningún notebook depende de esta carpeta.

## Por qué las quinasas

Son el caso de libro de una categoría **grande y promiscua**: 391 de los 14 083
dominios (2.8 %), pero tocan 13 921 de las 206 643 proteínas anotadas (6.7 %), y
están sobrerrepresentadas entre los blancos con bioactividad conocida. El
parámetro β del modelo penaliza exactamente eso —las categorías grandes—, así
que las quinasas son el grupo donde el efecto de β tiene que verse si existe.

## Cómo se produjo

Con `derivar.py` el 2026-09-21, sobre dos archivos privados:

| archivo | qué tiene | tamaño |
|---|---|---|
| `raw_data/targets/interpro.xml` | el catálogo de InterPro: 49 674 entradas con nombre y tipo | 246 MB |
| `acondicionarDB/mapeos/mapa_interpro.csv` | el diccionario de accession `IPR……` a código de `DB/` | 580 KB |

```bash
TDR_RAW=<carpeta raw_data> TDR_MAPEOS=<carpeta de mapeos> python derivar.py   # ~35 s
```

Se usa `interpro.xml` y no `protein2ipr.dat`, que trae la misma información
repetida por proteína y pesa 83 GB.

## Cómo se verificó

- El bloque que imprime `derivar.py` es **idéntico** al que está pegado en
  `comun/tdr.py`: mismas 391 entradas, mismo orden, mismos accessions.
- Los tres filtros se aplican en orden y cada uno se puede contar: 1 106 entradas
  de InterPro mencionan `kinase`, 412 de ellas son de tipo `Domain`, y 391
  existen en `DB/`. Las 21 que se caen son dominios que el recorte de la base no
  conserva.
- El tipo `Domain` es el mismo recorte que hace `tdr.cargar_db` (`IPR_DOMAIN`):
  marcar familias o superfamilias homólogas no cambiaría nada, porque no entran
  en la capa de anotaciones.
- `comun/tests/test_resultados.py` fija las 391 anotaciones y las 13 921
  proteínas, así que si la constante se toca sin querer, la prueba falla.

## Qué se revela y qué no

Se revela el accession de InterPro de **estos 391 dominios**, al lado de su
código. Es una decisión explícita (PROVENANCE.md §3.15) y es lo único del
diccionario privado que sale a la luz: los otros 36 533 accessions siguen sin
publicarse, y ningún mapeo de compuesto, blanco o cluster se toca.
