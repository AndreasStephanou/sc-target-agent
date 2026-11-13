#!/usr/bin/env bash
# Freeze the exact versions of the conda environment that produced a results tag.
#
#   bash scripts/freeze_env.sh [env-name]        # default: sctarget
#
# Writes environment/conda-linux-64.lock (conda packages by URL + md5; `conda create --file` installs
# them without a solver) and environment/pip-linux-64.lock.txt (the pip layer, `name==version`).
# Run it before every results tag so the numbers in docs/results_*.md are reproducible; environment.yml
# stays the loose, human-readable specification for a fresh install.
set -euo pipefail
ENV_NAME="${1:-sctarget}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="$ROOT/environment"
mkdir -p "$OUT"

PLATFORM="$(conda info --json | python3 -c 'import json,sys; print(json.load(sys.stdin)["platform"])')"
conda list -n "$ENV_NAME" --explicit --md5 > "$OUT/conda-$PLATFORM.lock"

TMP_JSON="$(mktemp)"
conda list -n "$ENV_NAME" --json > "$TMP_JSON"
python3 - "$TMP_JSON" "$OUT/pip-$PLATFORM.lock.txt" "$PLATFORM" <<'PY'
import datetime, json, sys
src, out, platform = sys.argv[1], sys.argv[2], sys.argv[3]
pk = json.load(open(src))
pypi = sorted((p["name"], p["version"]) for p in pk if p["channel"] == "pypi" and p["name"] != "sctarget")
lines = [
    "# pip layer of the sctarget environment, frozen by scripts/freeze_env.sh (conda list --json, channel pypi).",
    f"# Install after the conda layer (environment/conda-{platform}.lock), then `pip install --no-deps -e .`.",
    "# Every pip-installed package of the environment is listed, including tools unrelated to the pipeline.",
    f"# platform: {platform}; frozen: {datetime.date.today().isoformat()}",
]
lines += [f"{n}=={v}" for n, v in pypi]
open(out, "w").write("\n".join(lines) + "\n")
print(f"{len(pypi)} pip packages -> {out}")
PY
rm -f "$TMP_JSON"
echo "$(grep -c '^http' "$OUT/conda-$PLATFORM.lock") conda packages -> $OUT/conda-$PLATFORM.lock"
