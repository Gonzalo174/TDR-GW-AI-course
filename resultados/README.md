# resultados/

Salidas de los notebooks, una carpeta por análisis. Las crea `tdr.out(carpeta)`.

    resultados/<analisis>_out/            tablas NN_*.csv y el NN_meta.json de cada corrida
    resultados/<analisis>_out/figuras/    .png y .pdf

**Las tablas y los `NN_meta.json` se versionan** (53 MB): son la salida de la
corrida del 2026-09-10, y con ellas las figuras y las cifras de los entregables
se regeneran en segundos, sin volver a correr nada:

    make figuras        # las figuras, desde estas tablas (~40 s)
    make verificar      # pruebas + figuras + cifras del informe contra las tablas

**Las figuras no se versionan**: se dibujan desde estas tablas
(`comun/figuras.py`). Las que usa el informe se copian a `informe/figuras/`
con `make -C informe`, y esas sí se versionan.

Cada tabla lleva el número del notebook cuya celda de corrida la escribió
(CONVENCIONES.md §7). Volver a correr un notebook las reescribe.
