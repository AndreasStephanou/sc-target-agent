#!/usr/bin/env python3
"""Design check for a two-condition single-cell dataset, before any pipeline run.

Reads only `obs`, `var` and a slice of the count matrix (never the whole matrix), and prints:

  * where integer counts are (X, or which layer);
  * samples (libraries) per condition, cells per sample, and the imbalance between conditions;
  * pseudo-replication: if a donor/animal column is given, libraries per donor and whether the
    sample column is really the independent biological unit;
  * paired design: donors that appear in both conditions;
  * confounding: every sample-level categorical column (one value per sample) is cross-tabulated
    against the condition at the sample level; a column whose levels each occur in only one
    condition cannot be separated from the condition.

Nothing is modified. Exit code is 0 unless the file cannot be read; findings are printed as
WARNING lines and, with --json, written to a file so they can be quoted in a report.

Usage:
  python design_check.py data/my.h5ad --sample-key library --condition-key tissue \
      [--case Tumor] [--donor-key patient] [--max-levels 30] [--json design_check.json]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


def _read_obs(path: Path) -> pd.DataFrame:
    import anndata as ad  # obs/var only; X stays on disk

    a = ad.read_h5ad(path, backed="r")
    obs = a.obs.copy()
    layer_names = list(a.layers.keys()) if a.layers is not None else []
    try:
        a.file.close()
    except Exception:
        pass
    return obs, layer_names


def _looks_like_counts(h5: h5py.File, key: str, n_values: int = 200_000) -> dict:
    """Integrality check on a slice of X or a layer, without loading the matrix."""
    if key not in h5:
        return {"present": False}
    node = h5[key]
    if isinstance(node, h5py.Group):  # sparse: csr / csc
        data = node["data"]
        vals = data[: min(data.shape[0], n_values)]
        enc = node.attrs.get("encoding-type", b"")
        enc = enc.decode() if isinstance(enc, bytes) else str(enc)
    else:  # dense
        n_rows = min(node.shape[0], max(1, n_values // max(1, node.shape[1])))
        vals = np.asarray(node[:n_rows]).ravel()
        enc = "dense"
    vals = np.asarray(vals, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return {"present": True, "encoding": enc, "integer": None, "max": None}
    integer = bool(np.all(np.abs(vals - np.round(vals)) < 1e-6))
    return {
        "present": True,
        "encoding": enc,
        "integer": integer,
        "max": float(vals.max()),
        "nonneg": bool(vals.min() >= 0),
        "n_checked": int(vals.size),
    }


def _sample_level_columns(obs: pd.DataFrame, sample_key: str, max_levels: int) -> list[str]:
    """Categorical columns that take exactly one value per sample (library-level covariates)."""
    out = []
    for c in obs.columns:
        if c == sample_key:
            continue
        s = obs[c]
        if not (isinstance(s.dtype, pd.CategoricalDtype) or s.dtype == object or s.dtype == bool):
            continue
        if s.nunique(dropna=True) > max_levels or s.nunique(dropna=True) < 2:
            continue
        per_sample = obs.groupby(sample_key, observed=True)[c].nunique(dropna=True)
        if (per_sample <= 1).all():
            out.append(c)
    return out


def _cramers_v(ct: pd.DataFrame) -> float:
    n = ct.to_numpy().sum()
    if n == 0 or min(ct.shape) < 2:
        return float("nan")
    expected = np.outer(ct.sum(1), ct.sum(0)) / n
    with np.errstate(divide="ignore", invalid="ignore"):
        chi2 = np.nansum((ct.to_numpy() - expected) ** 2 / expected)
    k = min(ct.shape) - 1
    return float(np.sqrt(chi2 / (n * k))) if k > 0 else float("nan")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("h5ad")
    ap.add_argument("--sample-key", required=True, help="obs column: one value per library / biological replicate")
    ap.add_argument("--condition-key", required=True, help="obs column: condition")
    ap.add_argument("--case", default=None, help="condition value that is the case (disease / treated); others are reference")
    ap.add_argument("--donor-key", default=None, help="obs column naming the donor / animal, if libraries are nested in donors")
    ap.add_argument("--max-levels", type=int, default=30, help="ignore categorical columns with more levels than this")
    ap.add_argument("--json", default=None, help="write the findings to this JSON file")
    a = ap.parse_args(argv)

    path = Path(a.h5ad)
    if not path.exists():
        print(f"file not found: {path}", file=sys.stderr)
        return 2
    obs, layer_names = _read_obs(path)
    for k in (a.sample_key, a.condition_key):
        if k not in obs.columns:
            print(f"obs has no column {k!r}. Columns: {list(obs.columns)}", file=sys.stderr)
            return 2
    findings: dict = {"file": str(path), "n_cells": int(obs.shape[0]), "warnings": []}
    warn = findings["warnings"].append

    # ---- counts location
    with h5py.File(path, "r") as h5:
        counts = {"X": _looks_like_counts(h5, "X")}
        for ln in layer_names:
            counts[f"layers/{ln}"] = _looks_like_counts(h5, f"layers/{ln}")
    findings["counts"] = counts
    integer_sources = [k for k, v in counts.items() if v.get("present") and v.get("integer") and v.get("nonneg", True)]
    print(f"{obs.shape[0]} cells; matrices checked: " + "; ".join(
        f"{k}: {'integer' if v.get('integer') else 'non-integer'} (max {v.get('max')})" for k, v in counts.items() if v.get("present")))
    if not integer_sources:
        warn("no integer count matrix found in X or any layer: pseudobulk DESeq2 needs raw counts (set dataset.counts_layer, or rebuild from the raw deposit)")
    elif "X" not in integer_sources:
        print(f"  raw counts are in {integer_sources[0]} (set dataset.counts_layer accordingly); X is not integer")

    # ---- conditions and samples
    cond = obs[a.condition_key].astype(str)
    samp = obs[a.sample_key].astype(str)
    levels = sorted(cond.unique())
    findings["condition_levels"] = levels
    per_sample_cond = obs.groupby(a.sample_key, observed=True)[a.condition_key].nunique()
    mixed = per_sample_cond[per_sample_cond > 1]
    if len(mixed):
        warn(f"{len(mixed)} sample(s) carry more than one condition value ({list(mixed.index)[:5]}): the sample column is not a library/condition unit")
    sample_table = obs.groupby([a.condition_key, a.sample_key], observed=True).size()
    sample_table = sample_table[sample_table > 0]
    n_samples = sample_table.groupby(level=0, observed=True).size()
    cells_per_sample = sample_table.groupby(level=0, observed=True).agg(["min", "median", "max"])
    print(f"\ncondition levels: {levels}" + (f"  (case = {a.case})" if a.case else ""))
    for lv in levels:
        n = int(n_samples.get(lv, 0))
        cps = cells_per_sample.loc[lv] if lv in cells_per_sample.index else None
        print(f"  {lv}: {n} samples; cells per sample min/median/max = "
              + (f"{int(cps['min'])}/{int(cps['median'])}/{int(cps['max'])}" if cps is not None else "-"))
    findings["samples_per_condition"] = {k: int(v) for k, v in n_samples.items()}
    if len(levels) > 2:
        warn(f"{len(levels)} condition levels: the pipeline compares one case against the rest; set dataset.disease to one value and dataset.reference to the others, or split into pairwise runs")
    if a.case and a.case not in levels:
        warn(f"--case {a.case!r} is not a level of {a.condition_key}: {levels}")
    for lv, n in n_samples.items():
        if n < 2:
            warn(f"condition {lv!r} has {n} sample: replicate-aware tests are impossible; the pipeline should not be run on this design")
        elif n < 3:
            warn(f"condition {lv!r} has {n} samples: per-sample abundance and pseudobulk tests will be weak; every result must carry this limitation")
    if len(n_samples) >= 2:
        hi, lo = n_samples.max(), n_samples.min()
        if lo > 0 and hi / lo >= 3:
            warn(f"sample imbalance {int(hi)} vs {int(lo)}: the per-sample Mann-Whitney test is dominated by the larger group; report effect sizes alongside p-values")

    # ---- donors: pseudo-replication and pairing
    if a.donor_key:
        if a.donor_key not in obs.columns:
            warn(f"--donor-key {a.donor_key!r} not in obs")
        else:
            donor = obs[a.donor_key].astype(str)
            libs_per_donor = obs.groupby(a.donor_key, observed=True)[a.sample_key].nunique()
            donors_per_cond = obs.groupby(a.condition_key, observed=True)[a.donor_key].nunique()
            findings["donors_per_condition"] = {k: int(v) for k, v in donors_per_cond.items()}
            findings["libraries_per_donor_max"] = int(libs_per_donor.max())
            print(f"\ndonors per condition: {dict((k, int(v)) for k, v in donors_per_cond.items())}; libraries per donor max = {int(libs_per_donor.max())}")
            if (libs_per_donor > 1).any():
                multi = libs_per_donor[libs_per_donor > 1]
                shown = ", ".join(f"{k}: {int(v)}" for k, v in multi.head(5).items())
                warn(f"{len(multi)} donor(s) contribute more than one library ({shown}{', ...' if len(multi) > 5 else ''}): libraries from one donor are not independent replicates; the sample column over-counts replication. Either use the donor as sample_key (pool its libraries) or state the inflation in the report")
            both = set()
            for d, grp in obs.groupby(a.donor_key, observed=True):
                if grp[a.condition_key].nunique() > 1:
                    both.add(str(d))
            if both:
                findings["paired_donors"] = sorted(both)
                print(f"  paired design: {len(both)} donor(s) appear in both conditions ({sorted(both)[:6]})")
                print("  note: the within-cluster disease-vs-reference contrast is not blocked by donor; donor effects will inflate reference dispersion. State this in the write-up.")

    # ---- confounding: sample-level covariates vs condition
    cov_cols = _sample_level_columns(obs, a.sample_key, a.max_levels)
    cov_cols = [c for c in cov_cols if c != a.condition_key and c != a.donor_key]
    sample_meta = obs.groupby(a.sample_key, observed=True)[[a.condition_key] + cov_cols].first()
    findings["confounding"] = {}
    # drop aliases: a column with one level per sample carries no information beyond the sample
    # itself; a column whose partition of the samples equals the condition's is the condition renamed
    cond_part = sample_meta[a.condition_key].astype(str).factorize()[0]
    aliases = {}
    for c in list(cov_cols):
        col = sample_meta[c].astype(str)
        if col.nunique() == len(sample_meta):
            aliases[c] = "one level per sample (alias of the sample column)"
        elif col.nunique() == len(set(cond_part)) and pd.crosstab(col, cond_part).gt(0).sum(axis=1).eq(1).all():
            aliases[c] = "same partition as the condition column (alias of the condition)"
    cov_cols = [c for c in cov_cols if c not in aliases]
    findings["aliases"] = aliases
    if aliases:
        print("\ncolumns skipped as aliases: " + "; ".join(f"{k}: {v}" for k, v in aliases.items()))
    if cov_cols:
        print("\nsample-level covariates vs condition (one row per sample):")
    for c in cov_cols:
        ct = pd.crosstab(sample_meta[c].astype(str), sample_meta[a.condition_key].astype(str))
        v = _cramers_v(ct)
        # a level that occurs in only one condition
        exclusive = int((ct.gt(0).sum(axis=1) == 1).sum())
        perfectly = bool(exclusive == ct.shape[0] and ct.shape[0] >= 2)
        # with nearly one level per sample (breed, age, batch id) every level is trivially exclusive
        # and Cramer's V is inflated; report, do not warn
        sparse_levels = ct.shape[0] > max(2, len(sample_meta) // 2)
        findings["confounding"][c] = {"cramers_v": None if np.isnan(v) else round(v, 3), "levels": int(ct.shape[0]),
                                      "levels_exclusive_to_one_condition": exclusive, "perfectly_confounded": perfectly,
                                      "too_many_levels_for_n_samples": sparse_levels}
        if sparse_levels and not (perfectly and ct.shape[0] == 2):
            print(f"  {c}: {ct.shape[0]} levels for {len(sample_meta)} samples; association with condition cannot be estimated (V = {v:.2f} is inflated). Check by hand if it matters.")
            continue
        status = "PERFECTLY CONFOUNDED" if perfectly else ("strongly associated" if v >= 0.7 else "ok")
        print(f"  {c}: {ct.shape[0]} levels, Cramer's V = {v:.2f}, {exclusive} level(s) exclusive to one condition -> {status}")
        if perfectly:
            warn(f"covariate {c!r} is perfectly confounded with the condition (every level occurs in only one condition): condition and {c} cannot be separated; any result may be a {c} effect")
        elif v >= 0.7:
            warn(f"covariate {c!r} is strongly associated with the condition (Cramer's V {v:.2f}); check it before interpreting enrichment")
    if not cov_cols:
        print("\nno other sample-level categorical columns found to check for confounding (batch, chemistry, sex, site). If they exist outside the object, check them by hand.")

    print()
    for w in findings["warnings"]:
        print("WARNING:", w)
    if not findings["warnings"]:
        print("no design warnings")
    if a.json:
        Path(a.json).write_text(json.dumps(findings, indent=2, default=str))
        print(f"\nwrote {a.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
