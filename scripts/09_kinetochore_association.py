#!/usr/bin/env python3
"""Analysis C.4: kinetochore haplotypes beside the strain's centromere variants (descriptive).

Reads 07/08 outputs (kinetochore_genes.tsv), preprocessed.tsv (CDE calls) and
cdeii_steps.tsv (CDEII length vs S288C). S288C itself is the reference and is left out.

Outputs in config `output.dir`:
  kinetochore_cen_variants.tsv         strain x centromere: CDEI/CDEIII substitutions, CCG hit, CDEII delta
  kinetochore_strain_matrix.tsv        strain x (10 haplotypes + counts of centromeres differing from S288C)
  kinetochore_assoc_<gene>_<cde>.tsv   haplotype x centromere tables: CBF1-CDEI, CBF1-CDEII,
                                       CEP3-CDEIII, NDC10-CDEIII
  figures/kinetochore_haplotype_matrix.{pdf,png}
"""
import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_rgb
from matplotlib.patches import Patch, Rectangle

from centromere.cde import read_reference
from centromere.coevo import (cell_label, centromere_variants, haplotype_by_centromere,
                              reference_elements, strain_matrix)
from centromere.io import load_config, save_tsv

PAIRS = [("CBF1", "cdei"), ("CBF1", "cdeii"), ("CEP3", "cdeiii"), ("NDC10", "cdeiii")]
MAX_PLOT_STRAINS = 200
GENE_ORDER = ["CBF1", "CEP3", "NDC10", "CTF13", "SKP1", "CSE4", "SCM3", "MIF2", "OKP1", "AME1"]

# dataviz reference palette (light): categorical slots in fixed order, blue ramp, diverging blue/red
CATEGORICAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
NEUTRAL, SURFACE, INK, MUTED = "#f0efec", "#ffffff", "#0b0b0b", "#6f6e69"
CDE_ALT = "#c3c2b7"  # neutral, darker than NEUTRAL: CDEI differs (text gives the base)
BLUE = ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab"]       # 1, 2, 3, >=4 substitutions
RED = ["#f7d4d3", "#ef9a99", "#e34948", "#a92b2a"]        # +1-3, +4-7, +8-11, >=12 bp


def text_on(hex_color: str) -> str:
    r, g, b = to_rgb(hex_color)
    return INK if 0.299 * r + 0.587 * g + 0.114 * b > 0.6 else "#ffffff"


def hap_colors(genes: pd.DataFrame) -> dict[tuple[str, str], str]:
    """Shared (>= 2 strains) non-reference haplotypes get categorical slots in order of first
    appearance per gene; reference and singletons stay neutral (identity is in the text)."""
    out = {}
    for gene, g in genes.groupby("gene"):
        shared = [h for h, n in g["haplotype"].value_counts(sort=False).items() if h != "ref" and n >= 2]
        for k, h in enumerate(sorted(shared, key=lambda x: int(x[1:]))):
            out[(gene, h)] = CATEGORICAL[k % len(CATEGORICAL)]
    return out


def delta_color(d: float) -> str:
    if pd.isna(d) or d == 0:
        return NEUTRAL
    step = min(3, (abs(int(d)) - 1) // 4)
    return BLUE[step] if d < 0 else RED[step]


def plot(genes: pd.DataFrame, cen: pd.DataFrame, strains: list[str], order: list[str], path: Path,
         formats: list[str], dpi: int) -> None:
    colors = hap_colors(genes)
    g = genes.set_index(["strain", "gene"])
    # strains without a locus row for a centromere: a `?` cell (no CDE call)
    c = cen.set_index(["strain", "centromere"]).reindex(
        pd.MultiIndex.from_product([strains, order], names=["strain", "centromere"]))
    c["cde_ok"] = c["cde_ok"].fillna(False).astype(bool)
    panels = [("Kinetochore haplotypes", GENE_ORDER), ("CDEI vs S288C", order),
              ("CDEII length - S288C (bp)", order), ("CDEIII substitutions vs S288C", order)]
    widths = [len(cols) for _, cols in panels]
    fig, axes = plt.subplots(1, 4, figsize=(0.42 * sum(widths) + 2.5, 0.42 * len(strains) + 2.4),
                             gridspec_kw={"width_ratios": widths, "wspace": 0.06})
    for ax, (title, cols), k in zip(axes, panels, range(4)):
        for i, s in enumerate(strains):
            for j, col in enumerate(cols):
                if k == 0:
                    hap, status = g.loc[(s, col), "haplotype"], g.loc[(s, col), "status"]
                    fill = colors.get((col, hap), NEUTRAL)
                    txt = hap + ("\nfs" if status != "intact" else "")
                else:
                    row = c.loc[(s, col)]
                    if k == 1:
                        lab = cell_label(row, "cdei")
                        fill = NEUTRAL if lab in (".", "?") else CDE_ALT
                        txt = "" if lab == "." else lab.replace(",", "\n")
                    elif k == 2:
                        fill = delta_color(row["delta_s288c"])
                        txt = cell_label(row, "cdeii").replace(".", "")
                    else:
                        lab = cell_label(row, "cdeiii")
                        n = 0 if lab in (".", "?") else len(lab.split(","))
                        fill = BLUE[min(n, 4) - 1] if n else NEUTRAL
                        txt = "?" if lab == "?" else (str(n) if n else "")
                ax.add_patch(Rectangle((j + 0.04, i + 0.04), 0.92, 0.92, facecolor=fill, edgecolor=SURFACE, lw=1))
                if k == 3 and row["ccg_hit"] is True:
                    ax.add_patch(Rectangle((j + 0.1, i + 0.1), 0.8, 0.8, fill=False, edgecolor=INK, lw=1.6))
                if k == 0 and status != "intact":
                    ax.add_patch(Rectangle((j + 0.04, i + 0.04), 0.92, 0.92, fill=False, hatch="////",
                                           edgecolor=MUTED, lw=0))
                ax.text(j + 0.5, i + 0.5, txt, ha="center", va="center", fontsize=6.5 if "\n" in txt else 7.5,
                        color=text_on(fill))
        ax.set_xlim(0, len(cols))
        ax.set_ylim(len(strains), 0)
        ax.set_xticks(np.arange(len(cols)) + 0.5)
        ax.set_xticklabels([x.replace("CEN", "") if k else x for x in cols], rotation=90 if k == 0 else 0,
                           fontsize=8, color=INK)
        ax.set_yticks(np.arange(len(strains)) + 0.5)
        ax.set_yticklabels(strains if k == 0 else [], fontsize=8.5, color=INK)
        ax.tick_params(length=0)
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_title(title, fontsize=9.5, color=INK, loc="left")
        if k:
            ax.set_xlabel("centromere", fontsize=8, color=MUTED)

    handles = [Patch(facecolor=NEUTRAL, label="as S288C / unique haplotype"),
               Patch(facecolor=CATEGORICAL[0], label="haplotype shared by >= 2 strains (one colour per haplotype)"),
               Patch(facecolor=CDE_ALT, label="CDEI differs from S288C"),
               Patch(facecolor=NEUTRAL, hatch="////", edgecolor=MUTED, label="fs = frameshifted allele (collapsed assembly)"),
               Patch(facecolor=BLUE[2], label="CDEII shorter / CDEIII substitutions (darker = more)"),
               Patch(facecolor=RED[2], label="CDEII longer"),
               Patch(facecolor=SURFACE, edgecolor=INK, lw=1.6, label="CDEIII CCG triplet changed"),
               Patch(facecolor=NEUTRAL, label="? = no confident CDE call")]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=7.5, frameon=False,
               bbox_to_anchor=(0.5, -0.02))
    fig.subplots_adjust(left=0.05, right=0.995, top=0.88, bottom=0.2)
    path.parent.mkdir(parents=True, exist_ok=True)
    for fmt in formats:
        fig.savefig(path.with_suffix(f".{fmt}"), dpi=dpi, bbox_inches="tight", facecolor=SURFACE)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    config = load_config(args.config)
    out = Path(config["output"]["dir"])
    order = config["visualize"]["centromere_order"]

    ref = reference_elements(read_reference(config["reference"]))
    pre = pd.read_csv(out / "preprocessed.tsv", sep="\t")
    steps = pd.read_csv(out / "cdeii_steps.tsv", sep="\t")
    genes = pd.read_csv(out / "kinetochore_genes.tsv", sep="\t")
    genes = genes[(genes["strain"] != "S288C") & (genes["status"] != "absent")]

    cen = centromere_variants(pre, steps, ref)
    save_tsv(cen, out / "kinetochore_cen_variants.tsv")
    matrix = strain_matrix(cen, genes)[GENE_ORDER + ["n_cdei_alt", "n_cdei_cacatg", "n_cdeii_len_alt",
                                                     "n_cdeii_off_step", "n_cdeiii_alt", "n_ccg_hit",
                                                     "n_no_cde_call"]]
    save_tsv(matrix.reset_index(), out / "kinetochore_strain_matrix.tsv")

    pd.set_option("display.width", 250)
    for gene, el in PAIRS:
        t = haplotype_by_centromere(cen, genes, gene, el, order)
        save_tsv(t.reset_index(), out / f"kinetochore_assoc_{gene.lower()}_{el}.tsv")
        varied = [c for c in order if (t[c] != ".").any()]
        print(f"{gene} haplotype x {el.upper()} (centromeres with any difference from S288C):")
        print(t[[f"{gene}_haplotype", f"{gene}_status"] + varied].to_string())
        print()
    print(matrix.to_string())

    strains = sorted(matrix.index, key=lambda s: tuple(matrix.loc[s, ["CBF1", "CEP3", "NDC10"]]))
    vis = config["visualize"]
    if len(strains) > MAX_PLOT_STRAINS:   # one row per strain: unreadable, and > 2^16 px at 300 dpi
        print(f"{len(strains)} strains > {MAX_PLOT_STRAINS}: haplotype matrix figure skipped")
        return
    plot(genes, cen, strains, order, out / "figures" / "kinetochore_haplotype_matrix",
         vis["formats"], vis["dpi"])


if __name__ == "__main__":
    main()
