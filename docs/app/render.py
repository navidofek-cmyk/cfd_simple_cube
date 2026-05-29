"""Renderování Markdown obsahu na HTML (s podporou rovnic a zvýraznění kódu)."""
from __future__ import annotations
import markdown
from pathlib import Path

# Adresář s obsahem (../content vzhledem k tomuto souboru)
CONTENT_DIR = Path(__file__).resolve().parent.parent / "content"

MD_EXTENSIONS = [
    "extra",            # tabulky, fenced code, footnotes...
    "codehilite",       # zvýraznění syntaxe přes Pygments
    "toc",              # obsah + kotvy nadpisů
    "sane_lists",
    "admonition",       # poznámkové bloky (note/warning)
]
MD_EXT_CONFIG = {
    "codehilite": {"guess_lang": False, "css_class": "highlight"},
    "toc": {"permalink": True},
}


def render_markdown(md_text: str) -> tuple[str, str]:
    """Vrátí (html, toc_html). Rovnice $...$ / $$...$$ zpracuje MathJax v prohlížeči."""
    md = markdown.Markdown(extensions=MD_EXTENSIONS, extension_configs=MD_EXT_CONFIG)
    html = md.convert(md_text)
    toc = getattr(md, "toc", "")
    return html, toc


def load_page(rel_path: str) -> tuple[str, str, str]:
    """Načte .md soubor z content/ a vrátí (titulek, html, toc).

    rel_path např. 'theory/01-navier-stokes' (bez přípony).
    Titulek = první H1 nadpis (# ...), nebo název souboru.
    """
    md_file = (CONTENT_DIR / rel_path).with_suffix(".md")
    if not md_file.exists() or not md_file.resolve().is_relative_to(CONTENT_DIR):
        raise FileNotFoundError(rel_path)
    text = md_file.read_text(encoding="utf-8")
    title = rel_path.rsplit("/", 1)[-1]
    for line in text.splitlines():
        if line.startswith("# "):
            title = line[2:].strip()
            break
    html, toc = render_markdown(text)
    return title, html, toc
