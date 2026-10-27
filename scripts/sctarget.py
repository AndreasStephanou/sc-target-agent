"""sc-target-agent: one command for a new dataset.

    python scripts/sctarget.py inspect   data/my.h5ad --name mydata          # guess columns, write config + key template
    python scripts/sctarget.py calibrate --config config/mydata.yaml         # scripted model, no tokens; checks QC, states, decoys
    python scripts/sctarget.py run       --config config/mydata.yaml         # agent run (Claude Code login by default)
    python scripts/sctarget.py report    --ground-truth data/ground_truth/mydata.yaml runs/mydata_*   # results table

Each step prints what to look at before the next. `run` refuses to start if the config still has
scaffold placeholders, and warns (does not block) if no answer key is committed: the pipeline
produces hypotheses without one; the evaluator needs one.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PY = sys.executable


def _run(args: list[str]) -> int:
    print("$", " ".join(args))
    return subprocess.call(args, cwd=ROOT)


def _config_ok(path: str) -> list[str]:
    text = Path(path).read_text()
    bad = [tok for tok in ("SAMPLE_COLUMN", "CONDITION_COLUMN", "DISEASE", "CONTROL", "DESCRIBE") if tok in text]
    return bad


def cmd_inspect(a):
    return _run([PY, "scripts/new_dataset.py", a.h5ad, "--name", a.name] + (["--force"] if a.force else []))


def cmd_calibrate(a):
    bad = _config_ok(a.config)
    if bad:
        print(f"config still has placeholders {bad}: edit {a.config} first")
        return 2
    run_dir = a.run_dir or f"runs/{Path(a.config).stem}_calibrate"
    rc = _run([PY, "scripts/run_pipeline.py", "--config", a.config, "--model", "scripted", "--run-dir", run_dir])
    if rc:
        return rc
    print("\nCalibration run finished. Before an agent run, look at:")
    print(f"  - {run_dir}/report.md section 1 (QC: cells kept, per-type removal) and section 3 (which clusters were enriched, split candidates, doublet/stress suspects)")
    print("  - whether the states it recorded are the populations you expected, and whether anything looks like a dissociation artefact or a doublet cluster")
    if a.ground_truth:
        _run([PY, "scripts/run_eval.py", "--ground-truth", a.ground_truth, run_dir])
        print("  - the eval line: state found? decoys above the target? If the scripted floor already finds the answer, the agent run is testing what reasoning adds.")
    return 0


def cmd_run(a):
    bad = _config_ok(a.config)
    if bad:
        print(f"config still has placeholders {bad}: edit {a.config} first")
        return 2
    gt = Path("data/ground_truth") / f"{Path(a.config).stem}.yaml"
    if not gt.exists() or gt.stat().st_size < 50:
        print(f"note: no answer key at {gt}; the run will produce hypotheses but cannot be scored. If you have a known answer, write and COMMIT the key before this run.")
    run_dir = a.run_dir or f"runs/{Path(a.config).stem}_agent_1"
    args = [PY, "scripts/run_pipeline.py", "--config", a.config, "--model", a.model, "--run-dir", run_dir]
    if a.scoring_version:
        args += ["--scoring-version", str(a.scoring_version)]
    rc = _run(args)
    if rc == 0:
        print(f"\nDone. Read {run_dir}/report.md: section 5 (critic flags) first, then section 4 (ranking, with target_class and flags). Every number traces to {run_dir}/run.jsonl.")
        if gt.exists() and gt.stat().st_size >= 50:
            _run([PY, "scripts/run_eval.py", "--ground-truth", str(gt), run_dir])
    return rc


def cmd_report(a):
    return _run([PY, "scripts/collect_results.py", "--ground-truth", a.ground_truth, "--out", a.out] + a.run_dirs)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("inspect", help="inspect an .h5ad and scaffold config + key template"); p.add_argument("h5ad"); p.add_argument("--name", required=True); p.add_argument("--force", action="store_true"); p.set_defaults(fn=cmd_inspect)
    p = sub.add_parser("calibrate", help="scripted run (no tokens) + optional eval"); p.add_argument("--config", required=True); p.add_argument("--run-dir"); p.add_argument("--ground-truth"); p.set_defaults(fn=cmd_calibrate)
    p = sub.add_parser("run", help="agent run"); p.add_argument("--config", required=True); p.add_argument("--model", default="claude-code", help="claude-code (subscription login, default) or a Claude API model id (needs ANTHROPIC_API_KEY)"); p.add_argument("--run-dir"); p.add_argument("--scoring-version", type=int, choices=[1, 2], default=2); p.set_defaults(fn=cmd_run)
    p = sub.add_parser("report", help="results table from run directories"); p.add_argument("--ground-truth", required=True); p.add_argument("--out", default="docs/results.md"); p.add_argument("run_dirs", nargs="+"); p.set_defaults(fn=cmd_report)
    a = ap.parse_args(argv)
    sys.exit(a.fn(a))


if __name__ == "__main__":
    main()
