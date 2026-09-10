# resultados/

Salidas de los notebooks, una carpeta por análisis. Las crea `tdr.out(carpeta)`;
están acá versionadas vacías para que el repositorio se clone ya listo para correr.

    resultados/<analisis>_out/            tablas .csv y el NN_meta.json de cada corrida
    resultados/<analisis>_out/figuras/    .png y .pdf

`.gitignore` excluye el contenido: las tablas se regeneran corriendo las celdas
de corrida de los notebooks, y las figuras, a partir de las tablas y sin tocar
la base, con `python comun/figuras.py`. Las figuras que entran al informe se
copian a `informe/figuras/` (`make -C informe`), y esas sí se versionan.
