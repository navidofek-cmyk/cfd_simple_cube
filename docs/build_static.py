"""Generátor STATICKÉHO webu z dokumentace.

Vyrenderuje všechny Markdown stránky do samostatných HTML souborů (bez serveru),
přidá galerii výsledků s obrázky a vše uloží do docs/site/ — vhodné pro
GitHub Pages.

Spuštění:
    cd docs
    uv run python build_static.py
"""
from __future__ import annotations
import re
import shutil
from pathlib import Path
from jinja2 import Environment, FileSystemLoader

from app.render import load_page, render_markdown
from app.nav import NAV, flat_pages

BASE = Path(__file__).resolve().parent
SITE = BASE / "site"
IMG_SRC = BASE.parent          # kořen projektu — tam jsou *.png
TEMPLATES = BASE / "templates"
STATIC = BASE / "static"

# Obrázky výsledků (název souboru → popisek)
RESULT_IMAGES = [
    ("flow_visualization.png",   "Proudění — rychlost, tlak a proudnice (řezy XY a XZ)"),
    ("mesh_visualization.png",   "Výpočetní síť a poloha krychle v tunelu"),
    ("convergence.png",          "Historie konvergence (kontinuita + změna polí)"),
    ("matrix_sparsity.png",      "Struktura řídkých matic (sparsity pattern)"),
    ("matrix_coefficients.png",  "Prostorové rozložení koeficientů matic"),
    ("matrix_dense_block.png",   "Hustý blok matice kolem krychle (7bodový stencil)"),
]

# Statická navigace = běžná + sekce Výsledky
RESULTS_SLUG = "results"
STATIC_NAV = NAV + [("Výsledky", [(RESULTS_SLUG, "Vizualizace výsledků")])]


def slug_to_file(slug: str) -> str:
    """'theory/01-x' → 'theory__01-x.html', 'index' → 'index.html'."""
    if slug == "index":
        return "index.html"
    return slug.replace("/", "__") + ".html"


def rewrite_links(html: str) -> str:
    """Přepíše absolutní URL live aplikace na relativní statické soubory."""
    # /static/style.css → style.css, /static/menu.js → menu.js
    html = html.replace("/static/", "")
    # /page/<slug> → <flat>.html   (i s případnou kotvou #...)
    def repl_page(m):
        slug, anchor = m.group(1), m.group(2) or ""
        return f'href="{slug_to_file(slug)}{anchor}"'
    html = re.sub(r'href="/page/([^"#]+)(#[^"]*)?"', repl_page, html)
    # odkaz na domovskou "/" → index.html (jen přesně href="/")
    html = re.sub(r'href="/"', 'href="index.html"', html)
    return html


def build():
    # čistý výstupní adresář
    if SITE.exists():
        shutil.rmtree(SITE)
    SITE.mkdir(parents=True)
    (SITE / "img").mkdir()

    env = Environment(loader=FileSystemLoader(str(TEMPLATES)), autoescape=False)
    tpl = env.get_template("page.html")

    pages = [("index", "index")] + [(s, s) for s, _ in flat_pages() if s != "index"]
    slugs = [s for s, _ in pages]

    def prevnext(slug):
        all_slugs = [s for s, _ in flat_pages()] + [RESULTS_SLUG]
        if slug not in all_slugs:
            return None, None
        i = all_slugs.index(slug)
        labels = dict((s, l) for sec, items in STATIC_NAV for s, l in items)
        prev = (all_slugs[i-1], labels.get(all_slugs[i-1], "")) if i > 0 else None
        nxt  = (all_slugs[i+1], labels.get(all_slugs[i+1], "")) if i < len(all_slugs)-1 else None
        return prev, nxt

    # ── běžné stránky ────────────────────────────────────────────────
    count = 0
    for slug in slugs:
        title, content, toc = load_page(slug)
        prev, nxt = prevnext(slug)
        html = tpl.render(title=title, content=content, toc=toc,
                          nav=STATIC_NAV, active=slug, prev=prev, next=nxt)
        html = rewrite_links(html)
        (SITE / slug_to_file(slug)).write_text(html, encoding="utf-8")
        count += 1

    # ── stránka výsledků (galerie obrázků) ───────────────────────────
    md = ["# Vizualizace výsledků\n",
          "Výstupy řešiče (Re = 20). Obrázky generují Python skripty "
          "(`plot_fields.py`, `plot_mesh.py`, `plot_convergence.py`, "
          "`plot_matrix.py`) z dat řešiče.\n"]
    copied = 0
    for fname, caption in RESULT_IMAGES:
        src = IMG_SRC / fname
        if not src.exists():
            continue
        shutil.copy(src, SITE / "img" / fname)
        md.append(f"## {caption}\n")
        md.append(f"![{caption}](img/{fname})\n")
        copied += 1
    content, toc = render_markdown("\n".join(md))
    prev, nxt = prevnext(RESULTS_SLUG)
    html = tpl.render(title="Vizualizace výsledků", content=content, toc=toc,
                      nav=STATIC_NAV, active=RESULTS_SLUG, prev=prev, next=nxt)
    html = rewrite_links(html)
    (SITE / "results.html").write_text(html, encoding="utf-8")

    # ── statické soubory ─────────────────────────────────────────────
    shutil.copy(STATIC / "style.css", SITE / "style.css")
    shutil.copy(STATIC / "menu.js", SITE / "menu.js")
    # zabraň Jekyllu na GitHub Pages přepisovat soubory
    (SITE / ".nojekyll").write_text("")

    print(f"Hotovo: {count} stránek + galerie ({copied} obrázků) → {SITE}")
    print(f"Otevři: {SITE / 'index.html'}")


if __name__ == "__main__":
    build()
