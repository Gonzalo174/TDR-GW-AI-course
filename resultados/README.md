# resultados/

Salidas de los notebooks, una carpeta por análisis. Las crea `tdr.out(carpeta)`;
están acá versionadas vacías para que el repositorio se clone ya listo para correr.

    resultados/<analisis>_out/            tablas .csv y el NN_meta.json de cada corrida
    resultados/<analisis>_out/figuras/    .png y .pdf

`.gitignore` excluye el contenido: las salidas se regeneran corriendo los
notebooks. Las tablas y figuras que entran al informe se copian a `paper/`
(ver `paper/sincronizar.sh`), y esas sí se versionan.
