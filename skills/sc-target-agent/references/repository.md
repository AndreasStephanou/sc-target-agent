# The sc-target-agent repository

Code: <https://github.com/AndreasStephanou/sc-target-agent>. Tag `results-v3` is the version this
skill was written against (prompts v6, scoring v3, the lineage-relative attribution guard);
`results-v1` and `results-v2` predate the HCC and ovarian reviewer runs.

## Install and smoke test

```bash
git clone https://github.com/AndreasStephanou/sc-target-agent
cd sc-target-agent
conda env create -f environment.yml && conda activate sctarget
pytest -q                            # 87 passed, 1 skipped
python scripts/sctarget.py check-env # the sparse / large-data tool chain (what the tests and demo cannot exercise)
python scripts/sctarget.py demo      # offline synthetic demo: env check, then scripted PASS, baseline lets decoys through
```

## The commands, with their flags

```bash
python scripts/sctarget.py check-env
python scripts/sctarget.py demo
python scripts/sctarget.py inspect   data/my.h5ad --name mydata [--force]
python scripts/sctarget.py calibrate --config config/mydata.yaml [--ground-truth data/ground_truth/mydata.yaml] [--run-dir DIR] [--scoring-version 1|2|3]
python scripts/sctarget.py run       --config config/mydata.yaml [--model claude-code:<model id> | claude-code | <claude model id>] [--run-dir DIR] [--scoring-version 1|2|3]   # default scoring 3
python scripts/sctarget.py report    --ground-truth data/ground_truth/mydata.yaml [--out docs/results.md] runs/mydata_*

# data builders
python scripts/build_from_matrix.py --counts ... --meta ... --dry-run     # GEO 'one matrix + metadata' deposits -> .h5ad
python scripts/fetch_geo_10x.py GSEnnnnn --out data/x/raw --keep <title substring> ...   # GEO per-sample 10x deposits -> sample sheet (substring: 'tumor' matches 'peritumoral')
python scripts/build_from_10x_samples.py --sheet data/x/raw/samples.csv --out data/x/x.h5ad
python scripts/build_ccrcc.py                                              # worked example for per-sample 10x directories

# evaluation and re-scoring
python scripts/run_eval.py --ground-truth KEY [--k 5] [--stability] [--json out.json] runs/a runs/b runs/c   # --stability needs the key
python scripts/rescore.py --config CFG --scoring-version 1|2|3 [--recompute-attribution] [--full-candidates] [--suffix S] runs/...
python scripts/collect_results.py                                          # regenerates docs/results_*.md

# this skill's pre-run checks (in the skill's scripts/ directory, not the repository)
python -I <skill>/scripts/design_check.py data/my.h5ad --sample-key S --condition-key C [--case V] [--donor-key D] [--json out.json]
python -I <skill>/scripts/label_check.py  data/my.h5ad --sample-key S --condition-key C --genes G1 G2 [--layer counts] [--celltype-key K --celltype V]
```

`run` refuses to start while the config has scaffold placeholders and warns if no key is committed.
`--model claude-code:<model id>` pins the Claude Code CLI to a model; bare `claude-code` runs on
the CLI's own default, and in either case the run log records the model that answered
(`usage.model` on every model call, `facts.model_ids` in state.json; protocol 23r).
`--model <claude model id>` switches to the API backend (`ANTHROPIC_API_KEY` in the environment,
never in a file or chat). `SCTARGET_OFFLINE=1` keeps all lookups local (bundled fixtures; results
marked `source=fixture`).

## Config keys the skill relies on (`config/default.yaml` has every key with a comment)

| key | meaning |
|---|---|
| `dataset.sample_key`, `condition_key`, `celltype_key` | obs columns; `celltype_key: null` derives lineages from the data |
| `dataset.disease`, `dataset.reference` | the case value and the list of reference values (any two-condition design) |
| `dataset.condition_descriptions` | label -> plain-language meaning, shown to the model; write the design type here |
| `dataset.counts_layer` | layer holding raw counts (`null` if `X` is raw) |
| `dataset.exclude_samples`, `exclude_samples_reason` | pre-registered exclusion arm; both recorded in the run and the report |
| `dataset.forbidden_obs_columns` | identifiers the privacy guard must never let through |
| `dataset.species`, `human_symbol_var` | ortholog mapping for non-human data |
| `analysis.qc.*` | adaptive per-sample QC (`mode: adaptive`, `n_mads`, `max_mt_pct`, `mt_floor_pct`); tissue-specific caps |
| `analysis.lineages.ignore_annotation` | annotation-free arm |
| `analysis.attribution.min_reference_frac_in_cluster` (0.02), `min_reference_frac_of_condition` (0.005) | the thin-pool and lineage-representativeness guards |
| `analysis.differential_abundance.*`, `differential_expression.*` | enrichment and DE thresholds |
| `targets.objective` | `target` (default) or `marker` (response genes; druggability reported, not scored) |
| `targets.scoring_version`, `targets.weights`, `targets.exclude_gene_patterns`, `druggable_classes` | scoring |
| `privacy.*` | list and table caps, sample pseudonymisation |

## Tools and which role may call them

Analyst only: `run_qc`, `detect_doublets`, `integrate_and_cluster`, `subcluster`, `merge_clusters`,
`check_reference_libraries`. Nominator / critic: `map_orthologs`, `druggability_lookup`,
`literature_count`. All roles: `load_dataset`, `cluster_composition`, `differential_abundance`,
`find_markers`, `pseudobulk_de`, `derive_lineages`, `label_lineages`, `disease_programs`,
`condition_attribution`, `gene_specificity`, `enrichment`, `plot_genes`, `score_targets`,
`record_finding`. No role has a tool to exclude a library: that is a config decision.

## Where things are

| path | what |
|---|---|
| `sctarget/tools/` | the deterministic tools (QC, doublets, clustering, Harmony, lineages, differential abundance, pseudobulk DESeq2 with sibling / within-cluster / compartment / whole-lineage contrasts, specificity, attribution, druggability, scoring) |
| `sctarget/agents/` | roles and prompts (`PROMPT_VERSION` history at the top of `prompts.py`), orchestrator, model backends (API, scripted, Claude Code CLI), tool loop |
| `sctarget/privacy.py`, `runlog.py` | the privacy guard and the append-only audit log |
| `sctarget/eval/` | the evaluator (lenient/strict state recovery, target rank, decoys, novel, drug recall) and `stability.py` |
| `scripts/sctarget.py` | the CLI; the builders and evaluation scripts listed above |
| `config/` | per-dataset YAML; `default.yaml` documents every key |
| `data/ground_truth/` | frozen answer keys (never edited after a run); `caa.yaml`, `ccrcc.yaml`, `hcc.yaml`, `ovarian.yaml` are worked examples of key drafting |
| `docs/EVALUATION_PROTOCOL.md` | the method history: every rule and the named failure that motivated it (22 sections) |
| `docs/REPORT.md`, `docs/results_*.md` | the four-dataset write-up and the regenerated results tables |
| `skills/sc-target-agent/SKILL.md`, `.claude/skills/sc-target-agent/SKILL.md` | the repository's own copies of this skill (keep them in sync with the incubator copy) |
| `tests/` | 87 tests; must pass after any change |

## Versioning rules the skill depends on

- Tag `method-freeze-v1` is the analysis as frozen before the first model run.
- Every later change to prompts, scoring or analysis rules is a labelled version (`prompt_version`,
  `scoring_version`) recorded in each run's `config_used.json`, the report header and the eval line,
  and has a section in the protocol document naming the failure it fixed.
- Config-level changes (thresholds, caps, `exclude_samples` with a reason) are allowed per dataset
  and are recorded in the run; code-level changes are not made for one dataset.

## What leaves the machine

To Anthropic (Claude Code or API): the role prompts, tool descriptions, and tool results, which are
cluster-level tables (fold changes, enrichment statistics, marker lists), gene symbols and drug
names. The privacy guard rejects any payload containing per-cell data or real sample identifiers
(samples are pseudonymised before the model sees them). To DGIdb, Open Targets and NCBI
E-utilities: gene symbols. Nothing else; there are no uploads. The skill's helper scripts send nothing.
