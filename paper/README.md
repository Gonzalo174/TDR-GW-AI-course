# paper/ — informe LaTeX

Carpeta autocontenida: compila sin depender de nada fuera de `paper/`.

## Compilar

    make            # pdflatex -> bibtex -> pdflatex x2  => informe.pdf
    make clean      # borra auxiliares
    make sync       # re-copia figuras y tablas desde ../gon4/*_out

## Contenido

| Archivo/dir      | Qué es |
|------------------|--------|
| `informe.tex`    | Documento principal |
| `refs.bib`       | Bibliografía (BibTeX) |
| `figuras/`       | 28 figuras (`.pdf` vectorial + `.png`); LaTeX toma el `.pdf` |
| `datos/`         | Tablas CSV/JSON de las que salen todos los números del texto |
| `sincronizar.sh` | Re-copia `figuras/` y `datos/` desde `../gon4/*_out` |
| `Makefile`       | Reglas de compilación |

`datos/` reproduce la estructura de salida de los notebooks
(`analiceDB_out/`, `genome_prioritization_out/`, `huerfanas_out/`), incluidos
los `*_meta.json` con parámetros, versiones y fechas de corte de cada corrida.

**No copiadas** (>2 MB, no se citan valores individuales de ellas):
`analiceDB_out/02_tamanos_cluster.csv` (98 MB),
`analiceDB_out/03_promiscuidad_por_compuesto.csv` (5 MB),
`huerfanas_out/01_k1_drugs.csv` (2,5 MB).
Se regeneran con `make sync` si se bajara el umbral en `sincronizar.sh`.

## Dependencias LaTeX

Compila con TeX Live 2017 + `texlive-latex-extra` (lo que ya hay instalado).

Opcional, para guionado castellano correcto (hoy usa patrones ingleses y
rótulos definidos a mano, ver el bloque `\IfFileExists{spanish.ldf}` al
principio del `.tex`):

    sudo apt install texlive-lang-spanish

No hace falta ningún otro paquete: el documento evita `siunitx`, `biblatex`
y `xelatex` a propósito.
