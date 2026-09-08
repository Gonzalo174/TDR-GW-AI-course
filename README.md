# TDR-GW-AI-course

Proyecto final del curso GW-AI.

## Contenido

- `final-project.html` — consigna y página de presentación del proyecto.
- `DB/` — base de datos acondicionada del proyecto (CSV y CSV comprimidos).
  - `DB/meta.json` describe cada tabla: filas, columnas, tamaño, y el criterio
    de recorte aplicado al generarla.

## Datos

`DB/` fue generado por `acondicionar.py` a partir de la DB original de TDR.
El detalle completo (semilla, escala de peso, criterio de recorte, diccionarios
y esquema de cada tabla) está en `DB/meta.json`.

`01_edges_clusters_fingerprint` está partido en dos archivos `.gz`
(`part00`, `part01`) que deben concatenarse para reconstruir la tabla completa.
