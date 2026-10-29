"""Patch 37 (report only): replicate-level figures and the evidence table."""
from pathlib import Path

from sctarget.agents.orchestrator import run_pipeline


def test_report_has_library_level_figures_and_evidence_table(cfg, tmp_path):
    res = run_pipeline(cfg, model_name="scripted", run_dir=tmp_path / "rep")
    md = (res.run_dir / "report.md").read_text()
    names = sorted(p.name for p in (res.run_dir / "figures").glob("*.png"))
    assert any("library_expression" in n for n in names), names
    assert any("state_library_shares" in n for n in names), names
    assert not any("dotplot" in n or "umap_genes" in n for n in names), names
    assert "Evidence summary for the nominated candidates" in md
    assert "Marker-free library label check" in md
    assert "library_expression" in md and "state_library_shares" in md
    # every figure referenced in the markdown exists on disk
    import re
    for rel in re.findall(r"\]\((figures/[^)]+)\)", md):
        assert (res.run_dir / rel).exists(), rel
