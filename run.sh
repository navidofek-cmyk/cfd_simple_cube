#!/usr/bin/env bash
# Přeloží a spustí CFD řešič, volitelně i vizualizaci.
# Použití:
#   ./run.sh            # přeloží + spustí výpočet
#   ./run.sh plots      # po výpočtu spustí i Python vizualizace
#   ./run.sh docs       # spustí dokumentaci (FastAPI, port 8084)
set -euo pipefail

# Přepni se do adresáře projektu (kde leží tento skript)
cd "$(dirname "$(readlink -f "$0")")"

THREADS="${OMP_NUM_THREADS:-6}"

case "${1:-run}" in
  docs)
    exec ./docs/run.sh
    ;;
  plots)
    RUN_PLOTS=1
    ;;
  run)
    RUN_PLOTS=0
    ;;
  *)
    echo "Neznámý argument: $1  (použij: run | plots | docs)"; exit 1
    ;;
esac

echo "==> Překládám (make)…"
make

echo "==> Spouštím řešič na ${THREADS} vláknech…"
OMP_NUM_THREADS="${THREADS}" ./channelCube

if [ "${RUN_PLOTS}" = "1" ]; then
  echo "==> Vizualizace a validace…"
  python3 verify.py            || true
  python3 plot_fields.py       || true
  python3 plot_mesh.py         || true
  python3 plot_convergence.py  || true
  echo "==> Hotovo. Obrázky: flow_visualization.png, mesh_visualization.png, convergence.png"
fi
