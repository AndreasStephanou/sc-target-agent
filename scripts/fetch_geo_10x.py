"""Download per-sample 10x supplementary files for a GEO series and write a sample sheet.

    python scripts/fetch_geo_10x.py GSE270680 --out data/gastric/raw \
        --keep "Tumor" --keep "Normal" --drop "Lymph Node" --drop "Peripheral Blood"

For every GSM in the series whose title matches --keep (and none of --drop), lists the sample's
supplementary directory on the GEO FTP mirror, downloads the matrix / features (or genes) / barcodes
files, and appends a row to <out>/samples.csv with: sample (GSM), title, and the directory. The sheet
is what build_from_10x_samples.py reads; edit its condition / patient columns if the titles do not
parse (the script tries to pull them from the title).

Only GEO (ftp.ncbi.nlm.nih.gov) and NCBI E-utilities are contacted. Re-running skips files already present.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
import time
import urllib.request
from pathlib import Path

import requests

E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def gsm_list(gse: str) -> list[dict]:
    ids = requests.get(E + "esearch.fcgi", params={"db": "gds", "term": f"{gse}[ACCN] AND gsm[ETYP]", "retmax": 1000, "retmode": "json"}, timeout=60).json()["esearchresult"]["idlist"]
    out = []
    for i in range(0, len(ids), 60):
        s = requests.get(E + "esummary.fcgi", params={"db": "gds", "id": ",".join(ids[i:i + 60]), "retmode": "json"}, timeout=60).json()["result"]
        for k in ids[i:i + 60]:
            d = s.get(k, {})
            if str(d.get("accession", "")).startswith("GSM"):
                out.append({"gsm": d["accession"], "title": d.get("title", "")})
        time.sleep(0.35)
    return sorted(out, key=lambda r: r["gsm"])


def suppl_dir(gsm: str) -> str:
    stem = gsm[:-3] + "nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{stem}/{gsm}/suppl/"


def list_dir(url: str) -> list[str]:
    html = urllib.request.urlopen(url, timeout=60).read().decode("utf-8", "ignore")
    return sorted(set(re.findall(r'href="([^"]+\.gz)"', html)))


def download(url: str, dest: Path) -> None:
    if dest.exists() and dest.stat().st_size > 0:
        return
    tmp = dest.with_suffix(dest.suffix + ".part")
    urllib.request.urlretrieve(url, tmp)
    tmp.rename(dest)


def guess_condition(title: str, keep: list[str]) -> str:
    for k in keep:
        if k.lower() in title.lower():
            return k
    return ""


def guess_patient(title: str) -> str:
    m = re.search(r"(Donor|Patient|P)\s*_?(\d+)", title, re.I)
    return f"{m.group(1)}{m.group(2)}".replace(" ", "") if m else ""


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gse")
    ap.add_argument("--out", required=True)
    ap.add_argument("--keep", action="append", default=[], help="title must contain one of these as a substring, case-insensitive (repeatable). Mind substrings: 'tumor' also matches 'peritumoral'; check samples.csv before building")
    ap.add_argument("--drop", action="append", default=[], help="title must contain none of these (repeatable)")
    ap.add_argument("--dry-run", action="store_true", help="list what would be downloaded")
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for r in gsm_list(a.gse):
        t = r["title"]
        if a.keep and not any(k.lower() in t.lower() for k in a.keep):
            continue
        if any(d.lower() in t.lower() for d in a.drop):
            continue
        rows.append(r)
    print(f"{len(rows)} samples selected from {a.gse}")
    sheet = out / "samples.csv"
    with open(sheet, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sample", "condition", "patient", "title", "dir"])
        for r in rows:
            gsm, t = r["gsm"], r["title"]
            d = out / gsm
            cond, pat = guess_condition(t, a.keep), guess_patient(t)
            w.writerow([gsm, cond, pat, t, str(d)])
            try:
                files = list_dir(suppl_dir(gsm))
            except Exception as e:  # noqa: BLE001
                print(f"  {gsm}: could not list supplementary files ({e})", file=sys.stderr)
                continue
            wanted = [f for f in files if re.search(r"(matrix|features|genes|barcodes)", f, re.I)]
            print(f"  {gsm} [{cond or '?'} | {pat or '?'}] {t[:60]} -> {len(wanted)} files" + ("" if wanted else "  <-- NO 10x FILES; check the suppl dir"))
            if a.dry_run:
                continue
            d.mkdir(exist_ok=True)
            for f in wanted:
                download(suppl_dir(gsm) + f, d / f)
            time.sleep(0.3)
    print(f"\nwrote {sheet}. Check the condition / patient columns, then:\n  python scripts/build_from_10x_samples.py --sheet {sheet} --out data/<name>/<name>.h5ad")
    return 0


if __name__ == "__main__":
    sys.exit(main())
