#!/usr/bin/env python3
"""figuras.py — toda figura del repositorio, dibujada leyendo sólo `resultados/`.

Una corrida tiene dos mitades que cuestan muy distinto. Calcular las tablas
exige cargar la base (la capa química son 36 M de aristas) y a veces barrer una
grilla de parámetros: minutos. Dibujar una figura a partir de esas tablas lleva
segundos. Si las dos mitades están en la misma celda, retocar un eje obliga a
pagar la primera otra vez.

Por eso están separadas:

    notebook, celda de corrida   base -> tablas en resultados/<analisis>_out/
    este módulo                  tablas -> figuras en resultados/<analisis>_out/figuras/

Cada figura se declara con el decorador `@figura`, que registra de qué análisis
es y **qué tablas lee**. La función recibe un cargador que sólo le deja leer esas
tablas: si una figura intentara leer otra, o recalcular algo desde la base, falla
en vez de funcionar por accidente (lo vigila `comun/tests/test_figuras.py`). La
lista de tablas es además la procedencia de la figura: `--lista` la imprime.

    python comun/figuras.py                     # todas, ~30 s
    python comun/figuras.py huerfanas           # las de un análisis
    python comun/figuras.py 02_f01_roc_26       # una sola
    python comun/figuras.py --lista             # figura -> tablas que lee

Desde un notebook, la sección «Resultados» llama a `dibujar(nombre)`.
"""
from __future__ import annotations

import fnmatch
import importlib
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import pandas as pd

AQUI = Path(__file__).resolve().parent
if str(AQUI) not in sys.path:
    sys.path.insert(0, str(AQUI))
import tdr  # noqa: E402

# carpeta de análisis -> módulo con sus figuras
MODULOS = {
    "analiceDB": "figuras_analice",
    "genome_prioritization": "figuras_genome",
    "huerfanas": "figuras_huerfanas",
    "verificacion": "figuras_verificacion",
}


class FaltaTabla(FileNotFoundError):
    """La figura necesita una tabla que todavía no se produjo."""


@dataclass
class Figura:
    analisis: str
    nombre: str                  # NN_fMM_<nombre>, sin extensión
    lee: tuple                   # tablas (o patrones glob) de <analisis>_out/
    funcion: Callable
    doc: str = field(default="")

    @property
    def salidas(self) -> Path:
        return tdr.SALIDAS / f"{self.analisis}_out"

    @property
    def carpeta(self) -> Path:
        return self.salidas / "figuras"

    def faltan(self) -> list:
        """Tablas declaradas que no están en disco."""
        return [p for p in self.lee if not list(self.salidas.glob(p))]


REGISTRO: dict = {}


def figura(analisis: str, nombre: str, lee):
    """Registra una función de figura. `lee` son los nombres de las tablas de
    `resultados/<analisis>_out/` que la figura usa (se admiten patrones glob)."""
    if analisis not in MODULOS:
        raise ValueError(f"análisis desconocido: {analisis!r}")
    lee = (lee,) if isinstance(lee, str) else tuple(lee)

    def registrar(f):
        previa = REGISTRO.get(nombre)
        # la misma funcion vuelta a registrar es un `%autoreload`, no un duplicado
        if previa is not None and previa.funcion.__module__ != f.__module__:
            raise ValueError(f"figura duplicada: {nombre}")
        REGISTRO[nombre] = Figura(analisis, nombre, lee, f, (f.__doc__ or "").strip())
        return f
    return registrar


class Tablas:
    """Lo único que una figura puede leer: sus tablas declaradas."""

    def __init__(self, fig: Figura):
        self._fig = fig

    def _permitida(self, nombre):
        if not any(fnmatch.fnmatch(nombre, p) for p in self._fig.lee):
            raise KeyError(f"{self._fig.nombre} lee {nombre!r}, que no declaró en `lee` "
                           f"({', '.join(self._fig.lee)})")

    def ruta(self, nombre) -> Path:
        self._permitida(nombre)
        p = self._fig.salidas / nombre
        if not p.exists():
            raise FaltaTabla(f"falta {p.relative_to(tdr.RAIZ)}")
        return p

    def __call__(self, nombre, **kw) -> pd.DataFrame:
        """Una tabla. Vacía es un resultado, no un error: devuelve un DataFrame
        sin filas (pasa con las sugerencias de `huerfanas/04`)."""
        p = self.ruta(nombre)
        try:
            return pd.read_csv(p, **kw)
        except pd.errors.EmptyDataError:
            return pd.DataFrame()

    def varias(self, nb, patron, columna) -> pd.DataFrame:
        """Concatena `<nb>_<patron>.csv` con el sufijo como `columna` (ver
        `tdr.cargar_csvs`)."""
        if f"{nb}_{patron}.csv" not in self._fig.lee:     # el patrón mismo, declarado
            self._permitida(f"{nb}_{patron}.csv")
        return tdr.cargar_csvs(self._fig.salidas, nb, patron=patron, columna=columna)


# ---------------------------------------------------------------------------

def cargar_registro():
    """Importa los módulos de figuras de las cuatro carpetas (llenan REGISTRO)."""
    vacio = not REGISTRO          # p. ej. tras un `%autoreload` de este modulo
    for carpeta, modulo in MODULOS.items():
        d = str(tdr.RAIZ / carpeta)
        if d not in sys.path:
            sys.path.insert(0, d)
        if vacio and modulo in sys.modules:
            importlib.reload(sys.modules[modulo])
        else:
            importlib.import_module(modulo)
    return REGISTRO


def dibujar(nombre: str, guardar=True, plt=None):
    """Dibuja la figura `nombre` desde sus tablas y la guarda (pdf + png).

    Devuelve la figura, o None si no hay nada que dibujar (la función lo decide:
    una tabla de sugerencias vacía, por ejemplo)."""
    if not REGISTRO:
        cargar_registro()
    if nombre not in REGISTRO:
        raise KeyError(f"no hay figura registrada {nombre!r}")
    fig_def = REGISTRO[nombre]
    plt = plt or tdr.estilo(interactivo="matplotlib.pyplot" in sys.modules)
    fig = fig_def.funcion(Tablas(fig_def), plt)
    if fig is not None and guardar:
        tdr.guardar(fig, nombre, fig_def.carpeta)
    return fig


def regenerar(filtros=(), verbose=True) -> pd.DataFrame:
    """Dibuja todas las figuras que coinciden con `filtros` (análisis o nombre).

    Nunca carga la base. Las figuras cuyas tablas faltan se informan con el
    notebook que hay que correr, y no cuentan como error."""
    cargar_registro()
    plt = tdr.estilo(interactivo=False)
    filas = []
    for nombre, f in sorted(REGISTRO.items(), key=lambda kv: (kv[1].analisis, kv[0])):
        if filtros and not any(x in (f.analisis, nombre) for x in filtros):
            continue
        t0 = time.time()
        faltan = f.faltan()
        if faltan:
            estado, detalle = "faltan tablas", ", ".join(faltan)
        else:
            try:
                fig = dibujar(nombre, plt=plt) if verbose else _silencioso(nombre, plt)
                estado, detalle = ("ok", "") if fig is not None else ("sin datos", f.doc.split("\n")[0])
                if fig is not None:
                    plt.close(fig)
            except FaltaTabla as e:
                estado, detalle = "faltan tablas", str(e)
            except Exception as e:  # una figura rota no frena a las demás
                estado, detalle = "ERROR", f"{type(e).__name__}: {e}"
        filas.append({"analisis": f.analisis, "figura": nombre, "estado": estado,
                      "segundos": round(time.time() - t0, 1), "detalle": detalle})
    return pd.DataFrame(filas)


def _silencioso(nombre, plt):
    import contextlib
    import io
    with contextlib.redirect_stdout(io.StringIO()):
        return dibujar(nombre, plt=plt)


def lista() -> pd.DataFrame:
    cargar_registro()
    return pd.DataFrame([{"analisis": f.analisis, "figura": n, "lee": ", ".join(f.lee)}
                         for n, f in sorted(REGISTRO.items(),
                                            key=lambda kv: (kv[1].analisis, kv[0]))])


def main(args):
    pd.set_option("display.width", 200)
    pd.set_option("display.max_colwidth", 90)
    if "--lista" in args:
        print(lista().to_string(index=False))
        return 0
    t0 = time.time()
    r = regenerar(args, verbose=False)
    if r.empty:
        print(f"ninguna figura coincide con {args}; ver --lista")
        return 1
    print(r.to_string(index=False))
    print(f"\n{(r['estado'] == 'ok').sum()} de {len(r)} figuras en {time.time() - t0:.0f} s")
    if (r["estado"] == "faltan tablas").any():
        print("-> las que faltan se producen corriendo la celda de corrida del notebook "
              "cuyo número prefija la tabla")
    return 1 if (r["estado"] == "ERROR").any() else 0


if __name__ == "__main__":
    # Se delega en el módulo importado: los `figuras_*.py` registran en
    # `figuras.REGISTRO`, que no es el mismo objeto que el de `__main__`.
    import figuras
    sys.exit(figuras.main(sys.argv[1:]))
