"""FastAPI aplikace — servíruje dokumentaci CFD řešiče.

Spuštění:
    cd docs
    uv run uvicorn app.main:app --reload
"""
from __future__ import annotations
from pathlib import Path
from fastapi import FastAPI, Request, HTTPException
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates

from .render import load_page
from .nav import NAV, prev_next

BASE = Path(__file__).resolve().parent.parent
templates = Jinja2Templates(directory=str(BASE / "templates"))

app = FastAPI(title="CFD Cube — dokumentace", docs_url=None, redoc_url=None)
app.mount("/static", StaticFiles(directory=str(BASE / "static")), name="static")


def _render(request: Request, slug: str) -> HTMLResponse:
    try:
        title, html, toc = load_page(slug)
    except FileNotFoundError:
        raise HTTPException(404, f"Stránka '{slug}' neexistuje")
    prev, nxt = prev_next(slug)
    return templates.TemplateResponse(request, "page.html", {
        "title": title, "content": html, "toc": toc,
        "nav": NAV, "active": slug, "prev": prev, "next": nxt,
    })


@app.get("/", response_class=HTMLResponse)
def home(request: Request):
    return _render(request, "index")


@app.get("/page/{slug:path}", response_class=HTMLResponse)
def page(request: Request, slug: str):
    return _render(request, slug)
