# Compila el informe. `make` produce informe.pdf.
#
# Antes de compilar redibuja las figuras desde las tablas de resultados/ (sin
# tocar la base: `comun/figuras.py`, ~30 s), las trae y verifica que estén
# TODAS las que el .tex incluye: si falta una, LaTeX no
# falla, escribe el nombre del archivo en el lugar de la figura y el PDF sale
# con un hueco que se nota tarde.
TEX = informe

all: numeros figuras comprobar $(TEX).pdf presentacion.pdf

$(TEX).pdf: $(TEX).tex refs.bib numeros.tex tabla_optimos.tex
	pdflatex -interaction=nonstopmode $(TEX).tex >/dev/null
	-bibtex $(TEX) >/dev/null
	pdflatex -interaction=nonstopmode $(TEX).tex >/dev/null
	pdflatex -interaction=nonstopmode $(TEX).tex >/dev/null
	@echo "paginas: $$(pdfinfo $(TEX).pdf 2>/dev/null | awk '/Pages/{print $$2}' || echo '?')"
	@if grep -qE "Citation .* undefined|Reference .* undefined" $(TEX).log; then \
	  echo "AVISO: quedan citas o referencias sin resolver:"; \
	  grep -oE "(Citation|Reference) \`[^']*'" $(TEX).log | sort -u; \
	fi
	@if grep -q "??" numeros.tex; then \
	  echo "AVISO: hay numeros sin dato (??) — falta correr algun notebook:"; \
	  grep -o "newcommand{\\\\[A-Za-z]*}{??}" numeros.tex; \
	fi

# La presentacion de 5 minutos: mismas cifras (numeros.tex) y mismas figuras.
presentacion.pdf: presentacion.tex numeros.tex
	pdflatex -interaction=nonstopmode presentacion.tex >/dev/null
	pdflatex -interaction=nonstopmode presentacion.tex >/dev/null
	@if grep -qE "Reference .* undefined|Undefined control" presentacion.log; then \
	  echo "AVISO: presentacion con referencias o comandos sin resolver"; fi

numeros:
	@python ../informe/numeros.py

# Redibuja las figuras desde las tablas y las trae.
figuras:
	@python ../comun/figuras.py > /dev/null || python ../comun/figuras.py
	@mkdir -p figuras
	@for d in ../resultados/*_out/figuras; do cp -f $$d/*.pdf figuras/ 2>/dev/null || true; done
	@echo "figuras copiadas: $$(ls figuras/*.pdf 2>/dev/null | wc -l)"

# Toda figura que el .tex incluye tiene que existir.
comprobar:
	@falta=0; for f in $$(grep -ohE "figuras/[A-Za-z0-9_]+\.pdf" $(TEX).tex presentacion.tex | sort -u); do \
	  if [ ! -f "$$f" ]; then echo "FALTA la figura $$f"; falta=1; fi; done; \
	  if [ $$falta -eq 1 ]; then echo "-> correr los notebooks que las producen"; exit 1; fi; \
	  echo "todas las figuras del informe estan presentes"

limpiar:
	rm -f $(TEX).aux $(TEX).bbl $(TEX).blg $(TEX).log $(TEX).out $(TEX).toc \
	      presentacion.aux presentacion.log presentacion.nav presentacion.out presentacion.snm presentacion.toc

.PHONY: all numeros figuras comprobar limpiar
