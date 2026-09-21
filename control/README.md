# control/

`genoma_completo_v5/` son los óptimos por especie de la corrida v5, calculada de
forma independiente y con la métrica ya corregida. `genome_prioritization/`
compara sus resultados contra estos: es un control externo, no un insumo.

Un archivo por especie, nombrado por su **código** (`26.csv`): el código es la
clave con la que la especie entra en el modelo, y los nombres de archivo se
dejan como están para que la comparación contra la corrida v5 siga siendo byte a
byte. Qué organismo es cada código está en `tdr.SPECIES` (`26` es *Plasmodium
falciparum*). Cada archivo trae la grilla de parámetros con su `AUC01`, `AUC` y
corte de Youden.
