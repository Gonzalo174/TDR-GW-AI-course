# Compila el informe. Requiere pdflatex + bibtex (TeX Live).
TEX  = pdflatex -interaction=nonstopmode -halt-on-error
DOC  = informe

all: $(DOC).pdf

$(DOC).pdf: $(DOC).tex refs.bib
	$(TEX) $(DOC)
	-bibtex $(DOC)
	$(TEX) $(DOC)
	$(TEX) $(DOC)

sync:
	./sincronizar.sh

clean:
	rm -f $(DOC).aux $(DOC).log $(DOC).out $(DOC).toc \
	      $(DOC).bbl $(DOC).blg $(DOC).lof $(DOC).lot

distclean: clean
	rm -f $(DOC).pdf

.PHONY: all sync clean distclean
