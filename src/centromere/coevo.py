"""Analysis C.4: kinetochore haplotypes against the centromere variants of the same strain.

Descriptive only (n = 10): per strain x centromere, how CDEI / CDEII / CDEIII differ from
S288C, laid beside the strain's kinetochore haplotypes. No tests of association until
>= 50 genomes are available (spec C.4).
"""
import pandas as pd

from centromere.cde import CDEI_LEN, CDEIII_LEN

ELEMENTS = ("cdei", "cdeii", "cdeiii")


def substitutions(seq: str, ref: str) -> str:
    """`,`-joined 1-based substitutions of `seq` vs `ref` (equal length), e.g. `G6A`; `` if equal."""
    if len(seq) != len(ref):
        raise ValueError(f"length {len(seq)} != {len(ref)}")
    return ",".join(f"{r}{i}{s}" for i, (r, s) in enumerate(zip(ref, seq), 1) if r != s)


def reference_elements(reference: dict[str, str]) -> pd.DataFrame:
    """S288C CDEI (first 8 bp), CDEIII (last 26 bp), CDEII length and the CCG triplet of CDEIII
    (1-based start inside CDEIII; Cep3's Zn2Cys6 contacts it)."""
    rows = []
    for cen, s in reference.items():
        cdeiii = s[-CDEIII_LEN:]
        rows.append({"centromere": cen, "ref_cdei": s[:CDEI_LEN], "ref_cdeiii": cdeiii,
                     "ref_cdeii_length": len(s) - CDEI_LEN - CDEIII_LEN,
                     "ccg_start": cdeiii.find("CCG") + 1})
    return pd.DataFrame(rows)


def centromere_variants(pre: pd.DataFrame, steps: pd.DataFrame, ref: pd.DataFrame) -> pd.DataFrame:
    """Per strain x centromere: CDEI core (CACGTG / CACATG / ...), CDEI and CDEIII substitutions
    vs S288C, whether the CDEIII CCG is hit, CDEII length and its difference / helical-step
    class vs S288C. Rows with a failed CDE call (`cde_ok` False) keep NA for the CDEs."""
    d = pre[["strain", "centromere", "cdei_seq", "cdeiii_seq", "cde_ok"]].merge(ref, on="centromere")
    d["cdei_core"] = d["cdei_seq"].str[2:8]
    d["cdei_variants"] = [substitutions(s, r) for s, r in zip(d["cdei_seq"], d["ref_cdei"])]
    d["cdeiii_variants"] = [substitutions(s, r) for s, r in zip(d["cdeiii_seq"], d["ref_cdeiii"])]
    d["ccg_hit"] = [any(c != r for c, r in zip(s[p - 1:p + 2], q[p - 1:p + 2]))
                    for s, q, p in zip(d["cdeiii_seq"], d["ref_cdeiii"], d["ccg_start"])]
    for c in ("cdei_core", "cdei_variants", "cdeiii_variants", "ccg_hit"):
        d[c] = d[c].astype(object).where(d["cde_ok"])
    s = steps[["strain", "centromere", "cdeii_length", "delta_s288c", "class_s288c"]]
    d = d.merge(s, on=["strain", "centromere"], how="left")
    return d[["strain", "centromere", "cde_ok", "cdei_seq", "cdei_core", "cdei_variants",
              "cdeii_length", "delta_s288c", "class_s288c", "cdeiii_seq", "cdeiii_variants",
              "ccg_hit"]]


def cell_label(row: pd.Series, element: str) -> str:
    """Short text for one strain x centromere x element: `.` = as S288C, `?` = no CDE call."""
    if element == "cdeii":
        d = row["delta_s288c"]
        return "." if d == 0 else f"{int(d):+d}"
    if not row["cde_ok"]:
        return "?"
    return row[f"{element}_variants"] or "."


def haplotype_by_centromere(cen: pd.DataFrame, haplotypes: pd.DataFrame, gene: str,
                            element: str, order: list[str]) -> pd.DataFrame:
    """Rows = strains sorted by their `gene` haplotype, columns = centromeres; cells =
    `cell_label`. The first two columns are the haplotype and its allele status."""
    h = haplotypes[haplotypes["gene"] == gene].set_index("strain")[["haplotype", "status"]]
    lab = cen.assign(label=[cell_label(r, element) for _, r in cen.iterrows()])
    wide = lab.pivot(index="strain", columns="centromere", values="label")[order]
    out = h.join(wide, how="inner")
    key = out["haplotype"].map(lambda x: (x != "ref", x))  # ref first, then H1, H2...
    out = out.assign(_k=key).sort_values(["_k"]).drop(columns="_k")
    out.index.name = "strain"
    return out.rename(columns={"haplotype": f"{gene}_haplotype", "status": f"{gene}_status"})


def strain_matrix(cen: pd.DataFrame, haplotypes: pd.DataFrame) -> pd.DataFrame:
    """Strains x (haplotype of every kinetochore gene + per-element counts of centromeres
    differing from S288C)."""
    hap = haplotypes.pivot(index="strain", columns="gene", values="haplotype")
    ok = cen[cen["cde_ok"]]
    counts = pd.DataFrame({
        "n_cdei_alt": ok.groupby("strain")["cdei_variants"].apply(lambda x: int((x != "").sum())),
        "n_cdei_cacatg": ok.groupby("strain")["cdei_core"].apply(lambda x: int((x == "CACATG").sum())),
        "n_cdeii_len_alt": cen.groupby("strain")["delta_s288c"].apply(lambda x: int((x != 0).sum())),
        "n_cdeii_off_step": cen.groupby("strain")["class_s288c"].apply(lambda x: int((x == "off_step").sum())),
        "n_cdeiii_alt": ok.groupby("strain")["cdeiii_variants"].apply(lambda x: int((x != "").sum())),
        "n_ccg_hit": ok.groupby("strain")["ccg_hit"].apply(lambda x: int(x.astype(bool).sum())),
        "n_no_cde_call": cen.groupby("strain")["cde_ok"].apply(lambda x: int((~x.astype(bool)).sum())),
    })
    return hap.join(counts, how="inner")
