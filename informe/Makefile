# Compila el informe. `make` produce informe.pdf.
TEX = informe

all: $(TEX).pdf

$(TEX).pdf: $(TEX).tex refs.bib figuras
	pdflatex -interaction=nonstopmode $(TEX).tex >/dev/null
	-bibtex $(TEX) >/dev/null
	pdflatex -interaction=nonstopmode $(TEX).tex >/dev/null
	pdflatex -interaction=nonstopmode $(TEX).tex | tail -5

# Trae las figuras desde las salidas de los notebooks.
figuras:
	@mkdir -p figuras
	@for d in ../resultados/*_out/figuras; do cp -f $$d/*.pdf figuras/ 2>/dev/null || true; done
	@echo "figuras: $$(ls figuras/*.pdf 2>/dev/null | wc -l)"

limpiar:
	rm -f $(TEX).aux $(TEX).bbl $(TEX).blg $(TEX).log $(TEX).out $(TEX).toc

.PHONY: all figuras limpiar
