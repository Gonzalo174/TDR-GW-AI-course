#!/usr/bin/env bash
# Re-copia figuras y tablas desde las salidas de los notebooks (gon4/*_out)
# hacia esta carpeta, para que `paper/` sea autocontenido y compilable sola.
set -euo pipefail
cd "$(dirname "$0")"
SRC=../gon4
mkdir -p figuras datos
for d in "$SRC"/analiceDB_out "$SRC"/genome_prioritization_out "$SRC"/huerfanas_out; do
  cp -f "$d"/figuras/*.pdf figuras/ 2>/dev/null || true
  cp -f "$d"/figuras/*.png figuras/ 2>/dev/null || true
  p=$(basename "$d"); mkdir -p "datos/$p"
  # las tablas crudas >2 MB no se copian (ver README)
  find "$d" -maxdepth 1 -type f -size -2M -exec cp -f {} "datos/$p/" \;
done
echo "figuras: $(ls figuras | wc -l)  |  tablas: $(find datos -type f | wc -l)"
