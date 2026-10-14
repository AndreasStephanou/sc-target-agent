"""Score one or more runs against a ground-truth file; optionally report stability.

    python scripts/run_eval.py --ground-truth data/ground_truth/synthetic.yaml runs/run_a runs/run_b
    python scripts/run_eval.py --ground-truth data/ground_truth/synthetic.yaml --stability runs/rep_*
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from sctarget.eval import evaluate_run, format_eval, stability  # noqa: E402


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dirs", nargs="+")
    ap.add_argument("--ground-truth", required=True)
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--stability", action="store_true", help="also compute top-k / nomination stability across the given runs")
    ap.add_argument("--json", default=None, help="write all results to this json file")
    args = ap.parse_args(argv)

    results = []
    for d in args.run_dirs:
        res = evaluate_run(d, args.ground_truth, k=args.k)
        results.append(res.to_dict())
        print(format_eval(res))
        with open(Path(d) / "eval.json", "w") as fh:
            json.dump(res.to_dict(), fh, indent=2)
        print()
    out = {"runs": results}
    if args.stability and len(args.run_dirs) > 1:
        st = stability(args.run_dirs, k=10)
        out["stability"] = st
        print("Stability across runs")
        print("=" * 72)
        print(f"mean Jaccard top-10: {st['mean_jaccard_topk']}   mean Jaccard nominated: {st['mean_jaccard_nominated']}   rank-1 genes: {st['rank1_genes']}")
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=2)


if __name__ == "__main__":
    main()
