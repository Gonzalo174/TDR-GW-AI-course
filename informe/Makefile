# Compila el informe. `make` produce informe.pdf.
#
# Antes de compilar trae las figuras desde las salidas de los notebooks y
# verifica que estén TODAS las que el .tex incluye: si falta una, LaTeX no
# falla, escribe el nombre del archivo en el lugar de la figura y el PDF sale
# con un hueco que se nota tarde.
TEX = informe

all: numeros figuras comprobar $(TEX).pdf

$(TEX).pdf: $(TEX).tex refs.bib numeros.tex
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

numeros:
	@python ../informe/numeros.py

# Trae las figuras desde las salidas de los notebooks.
figuras:
	@mkdir -p figuras
	@for d in ../resultados/*_out/figuras; do cp -f $$d/*.pdf figuras/ 2>/dev/null || true; done
	@echo "figuras copiadas: $$(ls figuras/*.pdf 2>/dev/null | wc -l)"

# Toda figura que el .tex incluye tiene que existir.
comprobar:
	@falta=0; for f in $$(grep -oE "figuras/[A-Za-z0-9_]+\.pdf" $(TEX).tex | sort -u); do \
	  if [ ! -f "$$f" ]; then echo "FALTA la figura $$f"; falta=1; fi; done; \
	  if [ $$falta -eq 1 ]; then echo "-> correr los notebooks que las producen"; exit 1; fi; \
	  echo "todas las figuras del informe estan presentes"

limpiar:
	rm -f $(TEX).aux $(TEX).bbl $(TEX).blg $(TEX).log $(TEX).out $(TEX).toc

.PHONY: all numeros figuras comprobar limpiar
