# CFD Cube — dokumentace (FastAPI)

Kompletní dokumentace 3D CFD řešiče: teorie, C++ implementace, řešič.
Postavená na FastAPI, spravovaná přes `uv`.

## Spuštění

```bash
cd docs
uv run uvicorn app.main:app --host 0.0.0.0 --port 8084
```

Pak otevři <http://localhost:8084>.

Pro vývoj s automatickým restartem:

```bash
uv run uvicorn app.main:app --reload --port 8084
```

## Struktura

```
docs/
├── pyproject.toml      závislosti (uv)
├── app/
│   ├── main.py         FastAPI aplikace + routy
│   ├── render.py       Markdown → HTML (MathJax, zvýraznění)
│   └── nav.py          struktura navigace
├── content/            obsah v Markdownu
│   ├── theory/         9 stránek teorie
│   ├── cpp/            10 stránek C++ implementace
│   └── solver/         3 stránky řešiče
├── templates/          Jinja2 šablony
└── static/             CSS, JS
```

## Obsah

- **Teorie** — Navier-Stokes → FVM → SIMPLE → Rhie-Chow → CG → paralelizace
- **C++** — výklad každého souboru + „proč je tam co je" (OOP, šablony,
  operátory, RAII...)
- **Řešič** — SIMPLE smyčka, spuštění, validace

Obsah je v Markdownu — úpravy se projeví hned (s `--reload`).
