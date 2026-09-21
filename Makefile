# Puntos de entrada del repositorio. Tiempos medidos en una máquina de 48
# núcleos con 20 procesos (TDR_NCORE); ver README.md, "Tiempo de cómputo".
#
#   make verificar     pruebas + figuras desde las tablas + cifras del informe   ~1 min
#   make figuras       redibuja las figuras desde las tablas versionadas         ~40 s
#   make entregables   figuras, FIGURAS.md, informe, presentación y página       ~1,5 min
#   make corrida       corre los 10 notebooks en orden y reescribe las tablas    ~45 min
#
# `verificar` no escribe nada versionado: comprueba que las cifras que
# informe/numeros.py lee hoy de resultados/ sean las mismas que están en el
# informe versionado.

PY ?= python

NOTEBOOKS = analiceDB/01_dimensiones_red analiceDB/02_capas_y_conectividad \
            analiceDB/03_calidad_bioactividades \
            genome_prioritization/01_barrido_parametros \
            genome_prioritization/02_optimo_y_consistencia \
            genome_prioritization/03_comparacion_beta_promiscuidad \
            huerfanas/01_pseudohuerfanas huerfanas/02_cobertura_semilla \
            huerfanas/03_desorfanizacion_global huerfanas/04_aplicacion_especie \
            verificacion/10_equivalencia

.PHONY: verificar pruebas figuras entregables corrida

verificar: pruebas figuras
	@cp informe/numeros.tex /tmp/numeros_versionado.tex
	@$(PY) informe/numeros.py > /dev/null
	@if cmp -s informe/numeros.tex /tmp/numeros_versionado.tex; then \
	  echo "OK: las cifras del informe coinciden con las tablas de resultados/"; \
	else \
	  echo "DIFERENCIA entre las tablas y las cifras del informe:"; \
	  diff /tmp/numeros_versionado.tex informe/numeros.tex; \
	  cp /tmp/numeros_versionado.tex informe/numeros.tex; exit 1; \
	fi

pruebas:
	$(PY) -m unittest discover -s comun/tests -p "test_*.py"

figuras:
	$(PY) comun/figuras.py

entregables: figuras
	$(PY) comun/figuras.py --indice
	$(MAKE) -C informe
	$(PY) informe/generar_pagina.py

corrida:
	@for nb in $(NOTEBOOKS); do \
	  echo "== $$nb  $$(date +%H:%M:%S)"; \
	  (cd $$(dirname $$nb) && jupyter nbconvert --to notebook --execute --inplace \
	     --ExecutePreprocessor.timeout=-1 $$(basename $$nb).ipynb) || exit 1; \
	done
