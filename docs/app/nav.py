"""Struktura navigace — definuje stránky a jejich pořadí v menu."""

# Každá sekce: (titulek, [(slug, popisek), ...])
NAV = [
    ("Úvod", [
        ("index", "Přehled projektu"),
    ]),
    ("Teorie", [
        ("theory/01-navier-stokes",  "1. Navierovy–Stokesovy rovnice"),
        ("theory/02-fvm",            "2. Metoda konečných objemů"),
        ("theory/03-convection",     "3. Konvekce a upwind schéma"),
        ("theory/04-simple",         "4. Algoritmus SIMPLE"),
        ("theory/05-rhie-chow",      "5. Rhie–Chow interpolace"),
        ("theory/06-boundary",       "6. Okrajové podmínky"),
        ("theory/07-linear-system",  "7. Řídké matice a stencil"),
        ("theory/08-cg",             "8. Conjugate Gradient solver"),
        ("theory/09-parallel",       "9. Paralelizace (OpenMP)"),
    ]),
    ("C++ implementace", [
        ("cpp/00-overview",     "Architektura projektu"),
        ("cpp/01-vec3",         "Vec3.hpp — vektor a operátory"),
        ("cpp/02-field",        "Field.hpp — šablony a kontejnery"),
        ("cpp/03-mesh",         "Mesh.hpp — síť a maska"),
        ("cpp/04-fvmatrix",     "fvMatrix.hpp — matice a stencil"),
        ("cpp/05-cgsolver",     "CgSolver.hpp — lambda a workspace"),
        ("cpp/06-simple",       "Simple.cpp — jádro řešiče"),
        ("cpp/07-io",           "io.cpp — vstup/výstup"),
        ("cpp/08-main",         "main.cpp — vstupní bod"),
        ("cpp/09-language",     "Použité prvky C++ (reference)"),
    ]),
    ("Řešič", [
        ("solver/01-loop",       "SIMPLE smyčka krok po kroku"),
        ("solver/02-running",    "Spuštění a parametry"),
        ("solver/03-validation", "Validace a vizualizace"),
    ]),
]


def flat_pages():
    """Plochý seznam slugů v pořadí menu (pro prev/next navigaci)."""
    out = []
    for _, items in NAV:
        for slug, label in items:
            out.append((slug, label))
    return out


def prev_next(slug: str):
    pages = flat_pages()
    slugs = [s for s, _ in pages]
    if slug not in slugs:
        return None, None
    i = slugs.index(slug)
    prev = pages[i-1] if i > 0 else None
    nxt = pages[i+1] if i < len(pages)-1 else None
    return prev, nxt
