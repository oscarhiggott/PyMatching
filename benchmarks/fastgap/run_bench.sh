#!/usr/bin/env bash
# Regenerates every benchmark in this directory: DEMs, indices, per-shot CSVs, summaries and
# landmark tightness reports. Run from the repository root with the fastgap package installed.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
PY="${PYTHON:-python}"
CFG="$HERE/index_config.json"
cfg() { "$PY" -c "import json,sys; c=json.load(open('$CFG')); print(eval(sys.argv[1]))" "$1"; }

"$PY" "$HERE/gen_dems.py"
mkdir -p "$HERE/indices" "$HERE/results"
for name in $(cfg "' '.join(d['name'] for d in c['dems'])"); do
  dem="$HERE/dems/$name.dem"
  idx="$HERE/indices/$name.idx"
  "$PY" -m fastgap build-index --dem "$dem" --strategy "$(cfg "c['index']['strategy']")" \
      --num-landmarks "$(cfg "c['index']['num_landmarks']")" --out "$idx" > "$HERE/results/$name.index.json"
  "$PY" -m fastgap bench --dem "$dem" --idx "$idx" \
      --shots "$(cfg "c['bench']['shots']")" --warmup "$(cfg "c['bench']['warmup']")" \
      --threads $(cfg "' '.join(map(str, c['bench']['threads']))") --seed "$(cfg "c['bench']['seed']")" --no-upper-bound \
      --e2e --e2e-shots "$(cfg "c['bench']['e2e_shots']")" \
      --out "$HERE/results/$name.csv" --summary "$HERE/results/$name.summary.json" > /dev/null
  "$PY" -m fastgap landmark-report --dem "$dem" --idx "$idx" --shots "$(cfg "c['bench']['shots']")" \
      --seed "$(cfg "c['bench']['seed']")" --out "$HERE/results/$name.landmarks.json" > /dev/null
  echo "done $name"
done
