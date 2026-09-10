"""El entorno tiene todo lo que el repositorio importa.

Un paquete que falta no se nota hasta que alguien corre el notebook que lo usa
—`seaborn` y `networkx` aparecen en un solo import cada uno, lejos de la celda 1—
y para entonces la corrida ya se perdió. Este test recorre TODOS los `import` del
árbol, notebooks incluidos, y verifica que cada uno resuelva de verdad.
"""
import ast, json, pathlib, sys, importlib, unittest

RAIZ = pathlib.Path(__file__).resolve().parents[2]
LOCALES = {"tdr", "nucleo", "funciones_analice", "funciones_genome",
           "funciones_huerfanas", "funciones_verificacion",
           "figuras", "figuras_analice", "figuras_genome", "figuras_huerfanas",
           "figuras_verificacion", "numeros"}


def _sin_magics(src):
    """Las celdas con `%matplotlib` o `!pip` no parsean como Python."""
    return "\n".join("pass" if l.lstrip().startswith(("%", "!")) else l
                     for l in src.split("\n"))


def imports_del_repo():
    """{(modulo, nombre_importado)} de todos los .py y .ipynb del arbol."""
    encontrados, ilegibles = set(), []

    def leer(src, origen):
        try:
            arbol = ast.parse(_sin_magics(src))
        except SyntaxError as e:
            ilegibles.append((origen, str(e)))
            return
        for n in ast.walk(arbol):
            if isinstance(n, ast.Import):
                for a in n.names:
                    encontrados.add((a.name, ""))
            elif isinstance(n, ast.ImportFrom) and n.module and n.level == 0:
                for a in n.names:
                    encontrados.add((n.module, a.name))

    for p in sorted(RAIZ.rglob("*.py")):
        leer(p.read_text(), str(p))
    for p in sorted(RAIZ.rglob("*.ipynb")):
        for i, c in enumerate(json.loads(p.read_text())["cells"]):
            if c.get("cell_type") == "code":
                leer("".join(c["source"]), f"{p}[{i}]")
    return encontrados, ilegibles


class TestEntorno(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.imports, cls.ilegibles = imports_del_repo()
        for d in ("comun", "analiceDB", "genome_prioritization", "huerfanas", "verificacion"):
            if str(RAIZ / d) not in sys.path:
                sys.path.insert(0, str(RAIZ / d))

    def test_todas_las_celdas_parsean(self):
        """Una celda con sintaxis rota no falla hasta que se la ejecuta."""
        self.assertEqual(self.ilegibles, [],
                         f"código que no parsea: {[o for o, _ in self.ilegibles]}")

    def test_todo_lo_que_se_importa_esta_instalado(self):
        faltan = []
        for mod, nombre in sorted(self.imports):
            raiz = mod.split(".")[0]
            if raiz in sys.stdlib_module_names or raiz in LOCALES:
                continue
            try:
                m = importlib.import_module(mod)
                if nombre and not hasattr(m, nombre):
                    importlib.import_module(f"{mod}.{nombre}")
            except Exception as e:
                faltan.append(f"from {mod} import {nombre or '*'} ({type(e).__name__})")
        self.assertEqual(faltan, [], "falta instalar: " + "; ".join(faltan))

    def test_requirements_cubre_lo_que_se_usa(self):
        """Todo paquete externo importado tiene que estar fijado en requirements."""
        req = (RAIZ / "requirements.txt").read_text()
        fijados = {l.split("==")[0].strip().lower().replace("-", "_")
                   for l in req.split("\n") if "==" in l and not l.startswith("#")}
        alias = {"sklearn": "scikit_learn", "cv2": "opencv_python", "PIL": "pillow"}
        sin_fijar = set()
        for mod, _ in self.imports:
            raiz = mod.split(".")[0]
            if raiz in sys.stdlib_module_names or raiz in LOCALES:
                continue
            if alias.get(raiz, raiz).lower() not in fijados:
                sin_fijar.add(raiz)
        self.assertEqual(sin_fijar, set(), f"importados pero no fijados: {sin_fijar}")


if __name__ == "__main__":
    unittest.main()
