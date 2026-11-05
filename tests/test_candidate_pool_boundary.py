"""Protocol 23h: the scripted target role's candidate pool must keep a gene with one of the largest
significant effects even when it sits at significance rank 21-40, inside the stored 'stat' head but
outside the 20 the floor reads.

Motivating failure (ovarian floor, round 0.4b): FOLR1 was the 39th most significant gene of the
tumour cluster's sibling contrast (padj 5e-8) and the 9th largest effect (log2FC 4.8). The tool
stored it in the 40-gene 'stat' head, so it was not labelled 'effect'; the floor read the first 20
'stat' rows plus the 'effect' rows, and FOLR1 fell between the two.
"""
from __future__ import annotations

from sctarget.agents.model import ScriptedModel, ToolUse


class _Cfg:
    dataset = {"species": "homo_sapiens"}


class _Session:
    cfg = _Cfg()
    disease = "Tumor"

    def __init__(self, facts):
        self.facts = facts


def _rows():
    # 40 'stat' rows with modest effects (rank 1..40 by statistic), a top effect gene at stat rank 39,
    # then 20 'effect' rows, exactly as analysis._candidate_pool stores them
    rows = []
    for i in range(1, 41):
        g = "TARGET" if i == 39 else f"S{i:02d}"
        rows.append({"gene": g, "log2fc": 4.8 if g == "TARGET" else 1.2, "padj": 1e-9 * i, "pool": "stat"})
    for i in range(1, 21):
        rows.append({"gene": f"E{i:02d}", "log2fc": 4.0 - 0.1 * i, "padj": 1e-3, "pool": "effect"})
    return rows


def test_pool_keeps_large_effect_gene_inside_the_stat_head():
    facts = {
        "findings": {"disease_state": [{"cluster": "6", "basis": "composition", "genes": ["S01"]}]},
        "pseudobulk_de": {"6:cluster_vs_reference_clusters": {"top_up": _rows()}},
        "markers": {},
    }
    m = ScriptedModel(_Session(facts))
    m._facts = lambda: facts  # the scripted model reads the session facts through this hook
    steps = m._target_dynamic()
    spec = next(s for s in steps if isinstance(s, ToolUse) and s.name == "gene_specificity")
    genes = spec.input["genes"]
    assert "TARGET" in genes, "a gene with a top-10 effect at significance rank 39 must enter the pool"
    assert "S01" in genes and "S20" in genes, "the 20 most significant genes are kept"
    assert "E01" in genes, "the stored effect pool is kept"
    # the effect set is bounded: 20 head + at most 20 effects
    assert len(genes) <= 40
    # the weakest stored effect gene is displaced by the stronger in-head gene, not added on top
    assert "E20" not in genes
