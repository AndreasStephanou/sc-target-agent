#!/usr/bin/env bash
# Consolidation: apply the community package, regenerate every result from the run directories,
# tag, and (optionally) publish to GitHub for the CoPAIR community. Run block by block.
set -euo pipefail
cd /local/workdir/as3562/CoPAIR/sc-target-agent
conda activate sctarget 2>/dev/null || true

# 1. apply the community package
git apply fix29-community-package.patch
pytest -q tests/test_cli_and_collect.py
git add -A && git commit -m "community package: sctarget CLI, results collector, Claude Code skill, report draft"

# 2. regenerate the results tables from the run directories (no tokens)
python scripts/sctarget.py report --ground-truth data/ground_truth/caa.yaml   --out docs/results_caa.md   runs/caa_*/
python scripts/sctarget.py report --ground-truth data/ground_truth/ccrcc.yaml --out docs/results_ccrcc.md runs/ccrcc_*/
git add docs/results_caa.md docs/results_caa.json docs/results_ccrcc.md docs/results_ccrcc.json
git commit -m "results tables regenerated from run directories"
git tag results-v1

# 3. keep the repository small: run directories hold figures and tables; keep them out of git but back them up
echo "runs/" >> .gitignore; git add .gitignore; git commit -m "ignore runs/ (archived separately)" || true
tar czf ../sc-target-agent-runs-$(date +%Y%m%d).tgz runs/   # archive for the record (large; keep on the cluster, not in git)

# 4. publish (optional): Cornell GitHub or github.com. gh must be logged in.
#   gh repo create cornell-copair/sc-target-agent --private --source=. --remote=origin --push
#   or: gh repo create <your-user>/sc-target-agent --public --source=. --remote=origin --push
