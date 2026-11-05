"""Publication figure style for every figure the pipeline writes (tools and report alike).

Conventions follow the common life-science journal requirements (Nature-family style guides):
sans-serif type at 5-8 pt, single-column 89 mm / double-column 183 mm widths, no top or right
spines, thin axes, vector output with editable text (PDF font type 42, SVG text as text), and a
colour-blind-safe palette (Okabe & Ito 2008). Each figure is written three times with one stem:
`<stem>.png` (300 dpi, inlined in the HTML report), `<stem>.pdf` and `<stem>.svg`.

Nothing here changes any number; `apply_style()` only sets matplotlib rcParams.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

# Okabe & Ito (2008) qualitative palette
OKABE_ITO = {
    "orange": "#E69F00",
    "sky": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermilion": "#D55E00",
    "purple": "#CC79A7",
    "black": "#000000",
}
DISEASE_COLOR = OKABE_ITO["vermilion"]
REFERENCE_COLOR = OKABE_ITO["blue"]
NEUTRAL = "#8C8C8C"
LIGHT = "#D9D9D9"
ACCENT = OKABE_ITO["green"]
HIGHLIGHT = OKABE_ITO["orange"]
CATEGORICAL = [OKABE_ITO[k] for k in ("blue", "vermilion", "green", "orange", "purple", "sky", "yellow", "black")]
SEQUENTIAL_CMAP = "cividis"   # perceptually uniform, colour-blind safe
DIVERGING_CMAP = "RdBu_r"

MM = 1.0 / 25.4
SINGLE_COL = 89 * MM
MID_COL = 120 * MM
DOUBLE_COL = 183 * MM

_APPLIED = False


def apply_style() -> None:
    """Set matplotlib rcParams once per process."""
    global _APPLIED
    if _APPLIED:
        return
    import matplotlib

    matplotlib.use("Agg")
    from cycler import cycler

    matplotlib.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "Nimbus Sans", "DejaVu Sans"],
            "font.size": 7,
            "axes.titlesize": 8,
            "axes.titleweight": "bold",
            "axes.titlelocation": "left",
            "axes.labelsize": 7,
            "xtick.labelsize": 6,
            "ytick.labelsize": 6,
            "legend.fontsize": 6,
            "legend.title_fontsize": 6,
            "legend.frameon": False,
            "axes.linewidth": 0.6,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "lines.linewidth": 1.0,
            "lines.markersize": 3.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.grid": False,
            "axes.prop_cycle": cycler(color=CATEGORICAL),
            "figure.dpi": 110,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.02,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "image.cmap": SEQUENTIAL_CMAP,
        }
    )
    _APPLIED = True


def despine(ax: Any, left: bool = False, bottom: bool = False) -> None:
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if left:
        ax.spines["left"].set_visible(False)
    if bottom:
        ax.spines["bottom"].set_visible(False)


def panel_label(ax: Any, letter: str, x: float = -0.12, y: float = 1.05) -> None:
    """Bold lowercase panel letter outside the top-left corner, journal style."""
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="right")


def lineage_palette(names: list[str]) -> dict[str, str]:
    names = list(dict.fromkeys(str(n) for n in names))
    # lineage colours never reuse the condition colours (blue = reference, vermilion = disease)
    base = [OKABE_ITO[k] for k in ("green", "orange", "purple", "sky", "black")]  # no yellow: unreadable on white
    if len(names) <= len(base):
        return {n: base[i] for i, n in enumerate(names)}
    import matplotlib.pyplot as plt

    cmap = plt.get_cmap("tab20")
    return {n: (base[i] if i < len(base) else cmap((i - len(base)) % 20)) for i, n in enumerate(names)}


def save_all(fig: Any, stem: Path, formats: tuple[str, ...] = ("png", "pdf", "svg")) -> dict[str, str]:
    """Write the figure as PNG (300 dpi) + PDF + SVG with the same stem; returns {format: path}."""
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    out: dict[str, str] = {}
    for fmt in formats:
        p = stem.with_suffix(f".{fmt}")
        try:
            fig.savefig(p, format=fmt, dpi=300 if fmt == "png" else None)
            out[fmt] = str(p)
        except Exception:  # noqa: BLE001 - a vector backend failing must never lose the PNG
            if fmt == "png":
                raise
    return out
