# Para agentes

Si vas a corregir o reproducir este repositorio, empezá por
[`CORRECCION.md`](CORRECCION.md): tiene cada requisito de la consigna con su
evidencia y los valores que la verificación tiene que dar.

- Entorno: `conda create -n TDR_GW python=3.12 && pip install -r requirements.txt`
- Verificar sin correr nada (~2 min): `make verificar`
- Reproducir las tablas desde la base (~45 min, 20 procesos): `make corrida`
- No hace falta configurar rutas: todo se resuelve desde la raíz del repositorio.
- `analiceDB/03` necesita datos privados (`TDR_RAW`, `TDR_MAPEOS`); sin ellos
  `make corrida` lo saltea y se usa su tabla versionada.
- No escribas a mano en `informe/numeros.tex`, `informe/tabla_optimos.tex`,
  `FIGURAS.md` ni `index.html`: se generan desde las tablas.
