"""The community entry points: collect_results regenerates a table from run dirs; sctarget CLI guards."""
import subprocess
import sys
from pathlib import Path

from sctarget.agents.orchestrator import run_pipeline


def test_collect_results_from_run_dirs(cfg, ground_truth_path, tmp_path):
    r1 = run_pipeline(cfg, model_name="scripted", run_dir=tmp_path / "r1", write_report=False)
    out = tmp_path / "res.md"
    from scripts.collect_results import main as collect

    collect(["--ground-truth", str(ground_truth_path), "--out", str(out), str(r1.run_dir), str(tmp_path / "not_a_run")])
    md = out.read_text()
    assert "| r1 |" in md and "scripted" in md and "prompts v" in md
    assert out.with_suffix(".json").exists()


def test_cli_refuses_placeholder_config(tmp_path):
    cfg = tmp_path / "x.yaml"
    cfg.write_text("dataset:\n  path: x.h5ad\n  sample_key: SAMPLE_COLUMN\n  condition_key: c\n  disease: DISEASE\n  reference: [CONTROL]\n")
    root = Path(__file__).resolve().parents[1]
    r = subprocess.run([sys.executable, str(root / "scripts" / "sctarget.py"), "calibrate", "--config", str(cfg)], capture_output=True, text=True, cwd=root)
    assert r.returncode == 2 and "placeholders" in r.stdout
