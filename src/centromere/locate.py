"""L2 locator (S288C CEN flanks BLASTed to each assembly) and L3 locator (liftover of S288C
coordinates through minimap2 PAF alignments).

L2 (`locate_by_flanks`) works in BLAST coordinates: 1-based, inclusive.

Target = S288C R64 (the reference given to minimap2), query = assembly contig.
All coordinates here are 0-based, half-open, as in PAF.

The `cs` tag walks the target forward. On a reverse-strand block the query is
walked from its end: the k-th aligned query base is `qend - 1 - k`.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

_CS_OP = re.compile(r"(:\d+|\*[a-z][a-z]|[+\-][a-z]+|~[a-z]{2}\d+[a-z]{2})")
BASE_CODE = {"a": 1, "c": 2, "g": 3, "t": 4}
CODE_BASE = np.array(list("-ACGTN"))   # 0 = not covered, 5 = deleted / non-ACGT
DELETED = 5


def parse_cs(cs: str) -> pd.DataFrame:
    """Ops of a cs string with their start on target (t0) and query (q0), relative to the block.

    op: '=' match run, '*' substitution, '+' insertion in query, '-' deletion from query.
    """
    ops, tl, ql, base = [], [], [], []
    for tok in _CS_OP.findall(cs):
        c = tok[0]
        if c == ":":
            n = int(tok[1:])
            ops.append("="), tl.append(n), ql.append(n), base.append("")
        elif c == "*":
            ops.append("*"), tl.append(1), ql.append(1), base.append(tok[2])
        elif c == "+":
            ops.append("+"), tl.append(0), ql.append(len(tok) - 1), base.append(tok[1:])
        elif c == "-":
            ops.append("-"), tl.append(len(tok) - 1), ql.append(0), base.append(tok[1:])
        else:
            raise ValueError(f"unsupported cs op {tok!r} (spliced alignment)")
    if "".join(_CS_OP.findall(cs)) != cs:
        raise ValueError("cs string has unparsed characters")
    tl, ql = np.array(tl, dtype=np.int64), np.array(ql, dtype=np.int64)
    return pd.DataFrame({"op": ops, "tlen": tl, "qlen": ql, "base": base,
                         "t0": np.concatenate([[0], np.cumsum(tl)[:-1]]),
                         "q0": np.concatenate([[0], np.cumsum(ql)[:-1]])})


@dataclass
class Block:
    qname: str
    qlen: int
    qs: int
    qe: int
    strand: str
    tname: str
    ts: int
    te: int
    mapq: int
    cs: str = field(repr=False)
    _ops: pd.DataFrame | None = field(default=None, repr=False)

    @property
    def ops(self) -> pd.DataFrame:
        if self._ops is None:
            o = parse_cs(self.cs)
            if o["tlen"].sum() != self.te - self.ts or o["qlen"].sum() != self.qe - self.qs:
                raise ValueError(f"cs does not span block {self.qname}:{self.qs}-{self.qe}")
            self._ops = o
        return self._ops

    def qpos(self, qrel: int) -> int:
        """Query coordinate of the qrel-th aligned query base (in target orientation)."""
        return self.qs + qrel if self.strand == "+" else self.qe - 1 - qrel

    def lift(self, pos: int) -> int | str:
        """Query position of target position `pos`, or 'in_insertion_gap' if deleted in the query."""
        rel = pos - self.ts
        o = self.ops
        t = o[o["tlen"] > 0]
        i = int(np.searchsorted(t["t0"].to_numpy(), rel, side="right")) - 1
        r = t.iloc[i]
        if r["op"] == "-":
            return "in_insertion_gap"
        return self.qpos(int(r["q0"] + rel - r["t0"]))


def parse_paf(path, min_mapq: int = 20) -> list[Block]:
    """Primary (tp:A:P) blocks with MAPQ >= min_mapq and a cs tag."""
    blocks = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            tags = dict(t.split(":", 2)[::2] for t in f[12:])
            if tags.get("tp") != "P" or int(f[11]) < min_mapq or "cs" not in tags:
                continue
            blocks.append(Block(f[0], int(f[1]), int(f[2]), int(f[3]), f[4], f[5],
                                int(f[7]), int(f[8]), int(f[11]), tags["cs"]))
    return blocks


def _bridge(blocks: list[Block], chrom: str, pos: int, side: str, max_gap: int):
    """Query position for an unaligned target `pos` lying between two collinear blocks
    (same contig and strand, query in order, both gaps <= max_gap): the first query base
    after the left block (side='start') or the last one before the right block (side='end'),
    so a bridged interval contains the whole unaligned query segment."""
    left = [b for b in blocks if b.tname == chrom and b.te <= pos]
    right = [b for b in blocks if b.tname == chrom and b.ts > pos]
    if not left or not right:
        return "unaligned"
    L, R = max(left, key=lambda b: b.te), min(right, key=lambda b: b.ts)
    if L.qname != R.qname or L.strand != R.strand or R.ts - L.te > max_gap:
        return "unaligned"
    qgap = R.qs - L.qe if L.strand == "+" else L.qs - R.qe
    if qgap < 0 or qgap > max_gap:
        return "unaligned"
    if L.strand == "+":
        q = L.qe if side == "start" else R.qs - 1
    else:
        q = L.qs - 1 if side == "start" else R.qe
    return L.qname, q, L.strand


def liftover_paf(blocks: list[Block], chrom: str, pos: int, bridge: int = 0,
                 side: str = "start") -> tuple[str, int, str] | str:
    """(contig, query pos, strand) for target position `pos`, or the reason it does not lift:
    'unaligned', 'in_insertion_gap', 'multiple' (covered by more than one primary block).
    With bridge > 0 an unaligned position between collinear blocks is bridged (see _bridge)."""
    hits = [b for b in blocks if b.tname == chrom and b.ts <= pos < b.te]
    if not hits:
        return _bridge(blocks, chrom, pos, side, bridge) if bridge > 0 else "unaligned"
    if len(hits) > 1:
        return "multiple"
    b = hits[0]
    q = b.lift(pos)
    return q if isinstance(q, str) else (b.qname, q, b.strand)


def lift_interval(blocks: list[Block], chrom: str, start: int, end: int,
                  bridge: int = 0) -> tuple[str, int, int, str] | str:
    """Query interval (contig, qstart, qend, strand) for target [start, end), or a reason:
    the reason of a failed end, 'split_contigs', 'strand_flip' or 'order'."""
    a = liftover_paf(blocks, chrom, start, bridge, "start")
    b = liftover_paf(blocks, chrom, end - 1, bridge, "end")
    for x in (a, b):
        if isinstance(x, str):
            return x
    if a[0] != b[0]:
        return "split_contigs"
    if a[2] != b[2]:
        return "strand_flip"
    lo, hi = (a[1], b[1]) if a[2] == "+" else (b[1], a[1])
    if hi < lo:
        return "order"
    return a[0], lo, hi + 1, a[2]


def lift_window(blocks: list[Block], chrom: str, cen_start: int, cen_end: int, flank: int,
                step: int = 25, bridge: int = 0) -> dict:
    """Lift CEN +- flank; if an end does not lift onto the CEN's contig and strand, move it
    towards the CEN in `step`-bp steps until it does (largest lifted span containing the CEN)."""
    core = lift_interval(blocks, chrom, cen_start, cen_end, bridge)
    if isinstance(core, str):
        return {"status": core}
    contig, _, _, strand = core

    def ok(lo, hi):
        r = lift_interval(blocks, chrom, lo, hi, bridge)
        return not isinstance(r, str) and r[0] == contig and r[3] == strand

    left = flank
    while left > 0 and not ok(cen_start - left, cen_end):
        left = max(0, left - step)
    right = flank
    while right > 0 and not ok(cen_start - left, cen_end + right):
        right = max(0, right - step)
    _, qlo, qhi, _ = lift_interval(blocks, chrom, cen_start - left, cen_end + right, bridge)
    return {"status": "ok" if (left, right) == (flank, flank) else "shortened",
            "contig": contig, "qstart": qlo, "qend": qhi, "strand": strand,
            "left_kept": left, "right_kept": right}


def ref_alleles(blocks: list[Block], chrom: str, ref_seq: str) -> np.ndarray:
    """Allele code per target position of `chrom` (0 none, 1-4 ACGT, 5 deleted/non-ACGT),
    only where exactly one primary block covers the position; elsewhere 0."""
    n = len(ref_seq)
    cov = np.zeros(n, dtype=np.int16)
    out = np.zeros(n, dtype=np.int8)
    ref = np.frombuffer(ref_seq.upper().encode(), dtype=np.uint8)
    ref_code = np.full(256, DELETED, dtype=np.int8)
    for ch, c in BASE_CODE.items():
        ref_code[ord(ch.upper())] = c
    for b in blocks:
        if b.tname != chrom:
            continue
        cov[b.ts:b.te] += 1
        o = b.ops
        for op, tl, t0, base in zip(o["op"], o["tlen"], o["t0"], o["base"]):
            p = b.ts + t0
            if op == "=":
                out[p:p + tl] = ref_code[ref[p:p + tl]]
            elif op == "*":
                out[p] = BASE_CODE.get(base, DELETED)
            elif op == "-":
                out[p:p + tl] = DELETED
    out[cov != 1] = 0
    return out


FLANK_COLS = ["sample", "qseqid", "qstart", "qend", "qlen", "sseqid", "sstart", "send",
              "pident", "length", "evalue", "bitscore"]


def read_flank_hits(path) -> pd.DataFrame:
    """`blast.sh MODE=flanks` output with `centromere` and `side` (L/R) split from `CENn|L`."""
    h = pd.read_csv(path, sep="\t", header=None, names=FLANK_COLS)
    h[["centromere", "side"]] = h["qseqid"].str.split("|", expand=True)
    return h


def _flank_inner(h: pd.DataFrame) -> pd.DataFrame:
    """Strand and the CEN-facing end of each flank hit on the contig, extrapolated by the
    unaligned part of the flank (`ext`) so an HSP that stops short of the CEN still gives
    the CEN boundary. Flanks are in CEN orientation: L ends at CDEI, R starts after CDEIII."""
    h = h.copy()
    plus = h["send"] >= h["sstart"]
    left = h["side"] == "L"
    h["strand"] = np.where(plus, "+", "-")
    # L: inner = query end (qend -> send); R: inner = query start (qstart -> sstart)
    h["ext"] = np.where(left, h["qlen"] - h["qend"], h["qstart"] - 1)
    anchor = np.where(left, h["send"], h["sstart"])
    toward_cen = np.where(left == plus, 1, -1)    # L on + / R on -: CEN lies at higher coords
    h["inner"] = anchor + toward_cen * h["ext"]
    return h


def _refine_inner(anchor: pd.Series, other_inner: int, cands: pd.DataFrame, slack: int = 50) -> pd.Series:
    """The HSP of the same flank that lies closest to the CEN in the flank (smallest `ext`),
    on the anchor's contig and strand, collinear with it (its CEN-facing end not behind the
    anchor's by more than `slack` bp) and with its inner end still before the opposite flank.
    Recovers the boundary when an indel near the CEN splits the flank alignment."""
    toward = 1 if (anchor["side"] == "L") == (anchor["strand"] == "+") else -1
    a_end = anchor["inner"] - toward * anchor["ext"]           # aligned CEN-facing end of the anchor
    c = cands[(cands["sseqid"] == anchor["sseqid"]) & (cands["strand"] == anchor["strand"])
              & (cands["ext"] < anchor["ext"])]
    c_end = c["inner"] - toward * c["ext"]
    c = c[(toward * (c_end - a_end) >= -slack) & (toward * (other_inner - c["inner"]) > 0)]
    if c.empty:
        return anchor
    return c.sort_values(["ext", "bitscore"], ascending=[True, False]).iloc[0]


def locate_by_flanks(hits: pd.DataFrame, min_pident: float = 90, min_cov: float = 0.5,
                     max_locus: int = 20000, min_refine_len: int = 30) -> pd.DataFrame:
    """L2 locus per sample x CEN: the best left + right flank pair (by summed bitscore) on one
    contig, same strand, in CEN order, with the locus between their inner ends at most
    `max_locus` bp. Locus = 1-based inclusive [l2_start, l2_end]; `strand` = CEN orientation
    on the contig. Without a valid pair, `status` gives the reason:
    no_left / no_right / no_flanks, split_contigs, strand_flip, order (flanks overlap or are
    swapped) or too_far (> max_locus apart).

    The pair is chosen from HSPs covering >= `min_cov` of the flank; each inner end is then
    refined to the collinear HSP nearest the CEN (`_refine_inner`, any length >= `min_refine_len`),
    so an indel near the CEN does not shift the boundary. `ext_*` = flank bases between the
    used HSP and the CEN (extrapolated), `split_*` = the boundary came from a refining HSP."""
    allh = _flank_inner(hits[(hits["pident"] >= min_pident)
                             & (hits["qend"] - hits["qstart"] + 1 >= min_refine_len)])
    allh["hid"] = np.arange(len(allh))
    h = allh[(allh["qend"] - allh["qstart"] + 1) / allh["qlen"] >= min_cov]
    cand = dict(tuple(allh.groupby(["sample", "qseqid"])))
    rows = []
    keys = hits[["sample", "centromere"]].drop_duplicates()
    groups = dict(tuple(h.groupby(["sample", "centromere"])))
    for smp, cen in keys.itertuples(index=False):
        g = groups.get((smp, cen), h.iloc[:0])
        L, R = g[g["side"] == "L"], g[g["side"] == "R"]
        row = {"sample": smp, "centromere": cen, "n_left": len(L), "n_right": len(R)}
        if L.empty or R.empty:
            row["status"] = "no_flanks" if L.empty and R.empty else ("no_left" if L.empty else "no_right")
            rows.append(row)
            continue
        pairs = L.merge(R, on=["sseqid", "strand"], suffixes=("_l", "_r"))
        plus = pairs["strand"] == "+"
        pairs["start"] = np.where(plus, pairs["inner_l"], pairs["inner_r"]) + 1
        pairs["end"] = np.where(plus, pairs["inner_r"], pairs["inner_l"]) - 1
        pairs["len"] = pairs["end"] - pairs["start"] + 1
        valid = pairs[(pairs["len"] >= 0) & (pairs["len"] <= max_locus)]
        if valid.empty:
            bl, br = L.loc[L["bitscore"].idxmax()], R.loc[R["bitscore"].idxmax()]
            if bl["sseqid"] != br["sseqid"]:
                row["status"] = "split_contigs"
            elif bl["strand"] != br["strand"]:
                row["status"] = "strand_flip"
            else:
                row["status"] = "order" if pairs.loc[(pairs["sseqid"] == bl["sseqid"])
                                                     & (pairs["strand"] == bl["strand"]), "len"].min() < 0 \
                    else "too_far"
            rows.append(row)
            continue
        b = valid.loc[(valid["bitscore_l"] + valid["bitscore_r"]).idxmax()]
        al = allh[allh["hid"] == b["hid_l"]].iloc[0]
        ar = allh[allh["hid"] == b["hid_r"]].iloc[0]
        rl = _refine_inner(al, ar["inner"], cand[(smp, f"{cen}|L")])
        rr = _refine_inner(ar, rl["inner"], cand[(smp, f"{cen}|R")])
        lo, hi = (rl["inner"], rr["inner"]) if b["strand"] == "+" else (rr["inner"], rl["inner"])
        row.update(status="ok", contig=b["sseqid"], l2_start=int(lo + 1), l2_end=int(hi - 1),
                   l2_length=int(hi - lo - 1), strand=b["strand"],
                   pident_left=b["pident_l"], pident_right=b["pident_r"],
                   ext_left=int(rl["ext"]), ext_right=int(rr["ext"]),
                   split_left=bool(rl["hid"] != al["hid"]), split_right=bool(rr["hid"] != ar["hid"]))
        rows.append(row)
    return pd.DataFrame(rows)


def read_fasta(path) -> dict[str, str]:
    """Uppercase sequences by record id (.gz allowed)."""
    import gzip
    from pathlib import Path

    from Bio import SeqIO
    path = Path(path).expanduser()
    with (gzip.open(path, "rt") if path.suffix == ".gz" else open(path)) as fh:
        return {r.id: str(r.seq).upper() for r in SeqIO.parse(fh, "fasta")}


def assemblies(config: dict) -> pd.DataFrame:
    """accession, strain, assembly path and PAF path for every assembly in the sample map."""
    from pathlib import Path
    m = pd.read_csv(config["input"]["sample_map"], sep="\t")
    gdir = Path(config["cde"]["genome_dir"])
    pdir = Path(config["locate"]["paf_dir"])
    return m.assign(fasta=[gdir / f"{a}.fasta.gz" for a in m["accession"]],
                    paf=[pdir / f"{a}.paf" for a in m["accession"]])


def s288c_cens(config: dict) -> pd.DataFrame:
    """S288C CENs from the BED (0-based start), with the SGD sequence in CEN orientation."""
    from centromere.cde import read_reference
    bed = pd.read_csv(config["transposons"]["s288c_cen"], sep="\t", header=None,
                      names=["chrom", "start", "end", "centromere", "score", "strand"])
    ref = read_reference(config["reference"])
    return bed.assign(cen_seq=bed["centromere"].map(ref))


def compare_loci(df: pd.DataFrame, other: pd.DataFrame, tag: str, prefix: str,
                 ok: pd.Series | None = None) -> pd.DataFrame:
    """Join `other` (sample, centromere, contig, start, end) onto `df` (sample, centromere,
    contig, {prefix}_start, {prefix}_end) and add {tag}_contig/start/end, d_start_{tag},
    d_end_{tag} (df - other) and agree_{tag} (same contig, overlapping; NA where either side has
    no locus or `ok` is False)."""
    o = other.rename(columns={"contig": f"{tag}_contig", "start": f"{tag}_start", "end": f"{tag}_end"})
    df = df.assign(_ok=True if ok is None else ok.to_numpy(dtype=bool))   # positional, survives merge
    m = df.merge(o[["sample", "centromere", f"{tag}_contig", f"{tag}_start", f"{tag}_end"]],
                 on=["sample", "centromere"], how="left")
    m[f"d_start_{tag}"] = m[f"{prefix}_start"] - m[f"{tag}_start"]
    m[f"d_end_{tag}"] = m[f"{prefix}_end"] - m[f"{tag}_end"]
    ov = m[[f"{prefix}_end", f"{tag}_end"]].min(axis=1) - m[[f"{prefix}_start", f"{tag}_start"]].max(axis=1) + 1
    m[f"agree_{tag}"] = ((m["contig"] == m[f"{tag}_contig"]) & (ov > 0)).astype(object)
    bad = m[f"{tag}_contig"].isna() | m[f"{prefix}_start"].isna() | ~m["_ok"]
    m.loc[bad, f"agree_{tag}"] = pd.NA
    return m.drop(columns="_ok")


def l1_l2_l3_loci(config: dict) -> dict[str, pd.DataFrame]:
    """Existing loci per locator as sample, centromere, contig, start, end (1-based, inclusive);
    a locator whose output file is missing is left out."""
    from pathlib import Path
    out = Path(config["output"]["dir"])
    cols = ["sample", "centromere", "contig", "start", "end"]
    smap = pd.read_csv(config["input"]["sample_map"], sep="\t").set_index("strain")["accession"]
    res = {}
    if (out / "filtered.tsv").exists():
        f = pd.read_csv(out / "filtered.tsv", sep="\t")
        res["l1"] = f.rename(columns={"sseqid": "contig", "cen_start": "start", "cen_end": "end"})[cols]
    if (out / "l2_loci.tsv").exists():
        f = pd.read_csv(out / "l2_loci.tsv", sep="\t")
        f = f[f["status"] == "ok"].rename(columns={"l2_start": "start", "l2_end": "end"})
        res["l2"] = f[cols]
    if (out / "l3_loci.tsv").exists():
        f = pd.read_csv(out / "l3_loci.tsv", sep="\t")
        f = f[f["status"].isin(["lifted", "bridged"])].assign(sample=lambda d: d["strain"].map(smap))
        res["l3"] = f.rename(columns={"l3_start": "start", "l3_end": "end"})[cols]
    return res


# ---------------------------------------------------------------- phase 4: reconcile L1-L4

LOCATORS = ["l1", "l2", "l3", "l4"]
L4_FOUND = ["ok", "insertion", "cde", "length"]   # L4 statuses that place a locus (runner_up does not)
MAIN_TABLE = ["rescued", "insertion"]             # categories whose locus comes from L2/L3/L4


def locator_loci(config: dict, l1: pd.DataFrame) -> tuple[dict[str, pd.DataFrame], dict[str, set]]:
    """Loci of every locator as sample, centromere, contig, start, end (1-based, inclusive),
    strand (CEN orientation on the contig) and exact (False for an L3 locus bridged over an
    alignment gap, whose length is the gap's), plus the samples each locator ran on.

    `l1` is the main table (sseqid, cen_start, cen_end, strand); L2/L3/L4 are read from
    `l2_loci.tsv`, `l3_loci.tsv` and `<locate.hmm_dir>/hmm_loci.tsv`; a missing file = did not run."""
    from pathlib import Path
    out = Path(config["output"]["dir"])
    cols = ["sample", "centromere", "contig", "start", "end", "strand", "exact"]
    loci = {"l1": l1.rename(columns={"sseqid": "contig", "cen_start": "start", "cen_end": "end"})
            .assign(exact=True)[cols]}
    ran = {"l1": set(l1["sample"])}
    if (out / "l2_loci.tsv").exists():
        f = pd.read_csv(out / "l2_loci.tsv", sep="\t")
        ran["l2"] = set(f["sample"])
        loci["l2"] = f[f["status"] == "ok"].rename(columns={"l2_start": "start", "l2_end": "end"}) \
            .assign(exact=True)[cols]
    if (out / "l3_loci.tsv").exists():
        f = pd.read_csv(out / "l3_loci.tsv", sep="\t")
        smap = pd.read_csv(config["input"]["sample_map"], sep="\t").set_index("strain")["accession"]
        f["sample"] = f["strain"].map(smap)
        ran["l3"] = set(f["sample"].dropna())
        bed = pd.read_csv(config["transposons"]["s288c_cen"], sep="\t", header=None,
                          names=["chrom", "s", "e", "centromere", "score", "bed_strand"])
        f = f[f["status"].isin(["lifted", "bridged"])].merge(bed[["centromere", "bed_strand"]], on="centromere")
        # l3_strand is the alignment strand; a '-' CEN in S288C flips the CEN orientation
        f["strand"] = np.where((f["l3_strand"] == "+") == (f["bed_strand"] == "+"), "+", "-")
        loci["l3"] = f.rename(columns={"l3_start": "start", "l3_end": "end"}) \
            .assign(exact=f["status"] == "lifted")[cols]
    hmm = Path(config["locate"].get("hmm_dir", "")) / "hmm_loci.tsv"
    if config["locate"].get("hmm_dir") and hmm.exists():
        f = pd.read_csv(hmm, sep="\t")
        ran["l4"] = set(f["sample"])
        loci["l4"] = f[f["status"].isin(L4_FOUND)].rename(columns={"hmm_start": "start", "hmm_end": "end"}) \
            .assign(exact=True)[cols]
    return loci, ran


def _clusters(iv: pd.DataFrame) -> list[pd.DataFrame]:
    """Single-linkage groups of intervals that overlap on the same contig."""
    out = []
    for _, g in iv.sort_values(["contig", "start"]).groupby("contig", sort=False):
        cur, hi = [], None
        for i, r in g.iterrows():
            if cur and r["start"] > hi:
                out.append(g.loc[cur])
                cur = []
            hi = r["end"] if not cur else max(hi, r["end"])
            cur.append(i)
        out.append(g.loc[cur])
    return out


def _label(tags) -> str:
    return ",".join(t.upper() for t in LOCATORS if t in set(tags))


def reconcile(loci: dict[str, pd.DataFrame], ran: dict[str, set], ref_len: dict[str, int],
              insertion_bp: int = 30) -> pd.DataFrame:
    """One row per sample x CEN (every sample any locator ran on, every CEN in `ref_len`).

    Overlapping loci on one contig form a cluster; `support` = locators in it. Category:
      missing        no locator found the CEN
      conflict       loci on different contigs / not overlapping (listed in `clusters`, not resolved)
      insertion      an exact L2 or L3 locus is > `insertion_bp` longer than S288C
      rescued        L1 missed it, another locator found it
      sequence_only  L1 (and/or L4) found it, neither L2 nor L3 anchored it although one ran
      agree          every locator that ran on the sample found it
      partial        L1 and an anchor found it, some other locator did not (e.g. L4 runner_up)
    `source` = the locator whose coordinates the main table uses (L4 > L2 > L3 for rescued;
    the longest exact L2/L3 for insertion; L1 otherwise)."""
    tagged = pd.concat([d.assign(tag=t) for t, d in loci.items() if len(d)], ignore_index=True)
    by_key = dict(tuple(tagged.groupby(["sample", "centromere"])))
    samples = sorted(set().union(*ran.values()))
    rows = []
    for smp in samples:
        ran_here = [t for t in LOCATORS if smp in ran.get(t, set())]
        for cen, rl in ref_len.items():
            row = {"sample": smp, "centromere": cen, "ran": _label(ran_here)}
            iv = by_key.get((smp, cen))
            if iv is None:
                rows.append({**row, "category": "missing", "support": "", "n_clusters": 0})
                continue
            cl = _clusters(iv)
            row["n_clusters"] = len(cl)
            row["clusters"] = " | ".join(f"{_label(c['tag'])}:{c['contig'].iat[0]}:"
                                         f"{c['start'].min()}-{c['end'].max()}" for c in cl)
            for t in LOCATORS:
                x = iv[iv["tag"] == t]
                if len(x):
                    row[f"{t}_length"] = int((x["end"] - x["start"] + 1).max())
            if len(cl) > 1:
                l1c = [c for c in cl if "l1" in set(c["tag"])]
                c = l1c[0] if l1c else max(cl, key=len)
                row.update(category="conflict", support=_label(iv["tag"]))
            else:
                c = cl[0]
                sup = set(c["tag"])
                ex = c[c["tag"].isin(["l2", "l3"]) & c["exact"]]
                ins = ex[ex["end"] - ex["start"] + 1 > rl + insertion_bp]
                if len(ins):
                    cat = "insertion"
                elif "l1" not in sup:
                    cat = "rescued"
                elif not sup & {"l2", "l3"} and set(ran_here) & {"l2", "l3"}:
                    cat = "sequence_only"
                elif sup >= set(ran_here):
                    cat = "agree"
                else:
                    cat = "partial"
                row.update(category=cat, support=_label(sup))
                if cat == "insertion":
                    pick = ins.assign(n=ins["start"] - ins["end"]).sort_values(["n", "tag"]).iloc[0]
                    row["insertion_bp"] = int(pick["end"] - pick["start"] + 1 - rl)
                elif cat == "rescued":
                    pick = c.assign(o=c["tag"].map({"l4": 0, "l2": 1, "l3": 2})).sort_values("o").iloc[0]
                else:
                    pick = c[c["tag"] == "l1"].iloc[0]
                row.update(source=pick["tag"].upper(), contig=pick["contig"], start=int(pick["start"]),
                           end=int(pick["end"]), strand=pick["strand"])
                rows.append(row)
                continue
            p = c.sort_values("tag").iloc[0]
            row.update(source=p["tag"].upper(), contig=p["contig"], start=int(p["start"]),
                       end=int(p["end"]), strand=p["strand"])
            rows.append(row)
    df = pd.DataFrame(rows)
    front = ["sample", "centromere", "category", "support", "ran", "source", "contig", "start",
             "end", "strand", "n_clusters", "insertion_bp"]
    df = df.reindex(columns=front + [c for c in df.columns if c not in front])
    return df


def rescue_rows(rec: pd.DataFrame, genome_dir, reference: dict[str, str], strains: dict[str, str],
                flank: int = 150, radius: int = 20) -> pd.DataFrame:
    """Main-table rows (same CDE columns as `annotate_cde`) for `rescued` and `insertion` loci:
    CDEI/CDEIII placed with the S288C PWMs in the window given by the `source` locator."""
    from pathlib import Path

    from centromere.cde import CDEModel, _open_genome, locus_fields
    model = CDEModel(reference)
    todo = rec[rec["category"].isin(MAIN_TABLE)]
    rows = []
    for smp, g in todo.groupby("sample", sort=False):
        genome = _open_genome(Path(genome_dir), smp)
        for _, r in g.iterrows():
            f = locus_fields(genome[r["contig"]], int(r["start"]), int(r["end"]), r["strand"] == "+",
                             model, flank, radius)
            rows.append({"sample": smp, "centromere": r["centromere"], "sseqid": r["contig"],
                         "strain": strains.get(smp), "assembly": strains.get(smp), "blast_pieces": 0,
                         **f, "locator": r["source"]})
    return pd.DataFrame(rows)


def apply_reconcile(df: pd.DataFrame, rec: pd.DataFrame, rescued: pd.DataFrame) -> pd.DataFrame:
    """Main table after reconciliation: L1 rows of `insertion` loci are replaced by the rescued
    window, `rescued` rows are added; every row gets `locator` and `support`."""
    key = ["sample", "centromere"]
    ins = rec.loc[rec["category"] == "insertion", key]
    keep = df.merge(ins.assign(_drop=True), on=key, how="left")
    df = keep[keep["_drop"].isna()].drop(columns="_drop").assign(locator="L1")
    if len(rescued):
        df = pd.concat([df, rescued], ignore_index=True)
    return df.merge(rec[key + ["support", "category"]].rename(columns={"category": "locus_category"}),
                    on=key, how="left")
