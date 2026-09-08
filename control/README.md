# control/

`genoma_completo_v5/` son los óptimos por especie de la corrida v5, calculada de
forma independiente y con la métrica ya corregida. `genome_prioritization/`
compara sus resultados contra estos: es un control externo, no un insumo.

Un archivo por especie, nombrado por su **código** (`26.csv`), no por la especie:
la correspondencia código → organismo está en los diccionarios privados. Cada
archivo trae la grilla de parámetros con su `AUC01`, `AUC` y corte de Youden.
