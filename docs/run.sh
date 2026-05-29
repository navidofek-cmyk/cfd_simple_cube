#!/usr/bin/env bash
# Spustí dokumentační FastAPI server na portu 8084.
# Funguje odkudkoli — sám se přepne do složky docs/.
set -euo pipefail

PORT="${1:-8084}"

# Přepni se do adresáře, kde leží tento skript (docs/)
cd "$(dirname "$(readlink -f "$0")")"

echo "==> Synchronizuji závislosti (uv sync)…"
uv sync --quiet

# Uvolni port, pokud na něm něco visí
if command -v fuser >/dev/null 2>&1; then
    fuser -k "${PORT}/tcp" 2>/dev/null || true
    sleep 1
fi

IP=$(hostname -I 2>/dev/null | awk '{print $1}')
echo "==> Dokumentace poběží na:"
echo "      http://localhost:${PORT}"
[ -n "${IP}" ] && echo "      http://${IP}:${PORT}   (z jiného zařízení v síti)"
echo "    (zastavíš klávesami Ctrl+C)"
echo

exec uv run uvicorn app.main:app --host 0.0.0.0 --port "${PORT}"
