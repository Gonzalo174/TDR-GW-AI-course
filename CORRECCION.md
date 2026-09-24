# Guía para la corrección

Qué pide la consigna (`final-project.html`), dónde está la evidencia de cada
cosa y cómo comprobarla. Todo lo de esta página se puede verificar en unos
minutos, sin correr ningún modelo.

## Comprobarlo en dos minutos

```bash
conda create -n TDR_GW python=3.12 && conda activate TDR_GW
pip install -r requirements.txt
make verificar
```

`make verificar` hace tres cosas, y tiene que terminar diciendo
`OK: las cifras del informe coinciden con las tablas de resultados/`:

1. corre las 63 pruebas de `comun/tests/` (~25 s; 2 se saltean salvo con
   `TDR_TESTS_QUIMICA=1`, porque cargan la capa química entera);
2. redibuja las 36 figuras desde las tablas versionadas de `resultados/`, sin
   cargar la base (~40 s);
3. vuelve a calcular todas las cifras del informe desde esas tablas y comprueba
   que coincidan con las del informe versionado.

Valores que tiene que dar, y que las pruebas o el informe verifican:

| qué | valor esperado | dónde se comprueba |
|---|---|---|
| pruebas | 61 en verde, 2 salteadas | `make pruebas` |
| óptimos contra el control independiente | 16 de 16 iguales, ΔAUC01 = 0 | `test_resultados.test_los_optimos_reproducen_el_control` |
| comparación contra el oráculo v4 | 356 columnas, 0 discrepancias | `resultados/verificacion_out/10_equivalencia.csv` |
| figuras | 36 de 36 dibujadas | `make figuras` |

Reproducir las tablas desde la base, en cambio, lleva ~45 minutos con 20
procesos (`make corrida`; detalle en README, «Tiempo de cómputo»). Todos los
notebooks corren desde el repositorio. El único insumo calculado fuera de `DB/`
es la lista de compuestos promiscuos, que exige datos crudos privados: está
versionada en `datos_externos/promiscuidad/`, con su procedencia.

## La consigna, requisito por requisito

| la consigna pide | dónde está |
|---|---|
| **Un problema propio**, con muchos componentes | priorización de blancos terapéuticos, de la investigación del autor: datos (`DB/`), modelo (`comun/nucleo.py`), métrica corregida (`comun/tdr.py`), tres análisis, figuras, informe |
| **Hecho con agentes** | todo el repositorio se armó con Claude Code; los commits lo registran (`Co-Authored-By`). El camino, de v3 a v4 y de v4 a este repositorio, está en `informe/presentacion.pdf` (diapositivas 3 y 4), `PLAN.md` e `historia/` |
| **Procedencia de cada resultado**: qué archivo y qué función lo produjo, qué entró | cada tabla lleva el número del notebook que la escribió (CONVENCIONES §7) y su `NN_meta.json` (fecha, parámetros, fechas de las tablas de entrada, versiones); cada figura, en `FIGURAS.md`: notebook, función, tablas que lee |
| **Qué se escribió de cero y qué se llamó de una librería** | `PROVENANCE.md` §4 |
| **Cada elección con una alternativa defendible** | `PROVENANCE.md` §2 (datos) y §3 (port), cada fila con la alternativa descartada |
| **Cómo se verificó cada resultado**; «una figura con origen y sin verificación es sólo una afirmación» | `FIGURAS.md`, columna *verificación*, figura por figura, incluidas las que no tienen control independiente. Los tres niveles están en `PROVENANCE.md` §5: pruebas, control externo, oráculo |
| **Evidencia suficiente para defenderlo** | la comparación contra el oráculo v4 predice qué debe cambiar y qué no, y acierta (0 discrepancias); un error del modelo heredado de v4, corregido con una prueba (`PROVENANCE.md` §3.16); cuatro fallas silenciosas del port, encontradas por las pruebas o por la inspección de las salidas, en `PROVENANCE.md` §3 y en el informe §5 |
| **Una página HTML y un PDF** | `index.html`, `informe/informe.pdf` (y `informe/presentacion.pdf`) |
| **El repositorio, para máquinas: instrucciones, entorno, datos y verificaciones legibles por alguien que nunca habló con nosotros** | `README.md` (instalar, verificar, reproducir, tiempos); `requirements.txt` (versiones exactas, vigilado por `test_entorno.py`); `DB/` (la base, con `DB/meta.json`); `comun/tests/`; `Makefile` |

## Lo que no está, dicho explícitamente

- **La base original y su acondicionamiento**: `DB/` se deriva de datos que no
  se publican. Lo que se publica es la salida y el registro de cómo se produjo
  (`DB/meta.json`, `PROVENANCE.md` §2).
- **El cálculo de la lista de compuestos promiscuos**: necesita datos crudos
  privados; se versiona su salida (`datos_externos/promiscuidad/README.md`).
- **Los nombres de las anotaciones y la identidad de los organismos**: por
  diseño. Tres de los 16 organismos son identificables por su tipo, y está
  anotado.
- **Un control independiente del modelo**: el control de v5 corre el mismo
  `nucleo.py`, así que valida el port y no el modelo (CONVENCIONES §3).
- **La página todavía no está publicada**: se espera el contacto de los
  organizadores del curso.

Los pendientes están en `PLAN.md`, «Lo que queda».

## Para leer el código

| archivo | qué es |
|---|---|
| `CONVENCIONES.md` | las reglas de organización que cita el código (`CONVENCIONES.md §N`) |
| `comun/tdr.py` | rutas, carga de la base, métricas, estilo |
| `comun/nucleo.py` | el modelo; congelado salvo dos cambios, cada uno con su registro (PROVENANCE §3.4 y §3.16) |
| `comun/figuras.py` | el registro de figuras y su regeneración |
| `<analisis>/funciones_*.py` | el cálculo de cada análisis |
| `<analisis>/figuras_*.py` | las figuras de cada análisis |
| `<analisis>/NN_*.ipynb` | los notebooks, ya ejecutados, con sus salidas |
