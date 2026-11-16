#!/usr/bin/env python3
"""Download genome assemblies and reference files for the centromere pipeline.

    python scripts/download_data.py --set reference
    python scripts/download_data.py --set scrap142    [--limit N | --strains A,B] [--archive PATH]
    python scripts/download_data.py --set loegler1086 [--limit N | --strains A,B] [--archive PATH]

reference    S288C R64 genome (UCSC sacCer3, chrI..chrM) -> data/external/s288c/sacCer3.fa,
             SGD GFF -> data/external/sgd/, Peter 2018 and Loegler 2025 supplementary tables -> data/,
             Helsen 2026 Supplementary Data 4 and 6 -> data/external/helsen2026/.
scrap142     ScRAP, O'Donnell et al. 2023 (Zenodo 8210254, nuclear.tar.gz, 1.7 GB). One nuclear
             assembly per strain: `collapsed` where it exists, else `HP0` (142 strains).
loegler1086  Loegler et al. 2025 (Zenodo 15698884, GenomeAssemblies.tar.gz, 16.8 GB). The final
             assembly of each isolate as listed in Table S2 (`Filename`); phased isolates: HP1
             (1086 isolates). SecondBest assemblies and annotations are skipped.

`--haplotypes all` also keeps HP1+HP2 (ScRAP) / HP2 (Loegler) as extra samples; they double-count
diploids in population statistics, so the default is one assembly per strain.

The archive is streamed and never stored (use `--archive` for a local copy). Each kept assembly is
written gzipped to data/assemblies_<set>/<sample>.fasta.gz (via .tmp + rename; existing files are
skipped, so an interrupted download can be restarted). When the whole archive was read, its md5 is
checked against Zenodo's. The sample table data/samples_<set>.tsv (accession, assembly, strain;
accession = file stem) is what `input.sample_map` in config.<set>.yaml reads.
"""
import argparse
import gzip
import hashlib
import json
import re
import shutil
import sys
import tarfile
import urllib.request
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent

ZENODO = {
    "scrap142": ("8210254", "nuclear.tar.gz"),
    "loegler1086": ("15698884", "GenomeAssemblies.tar.gz"),
}
REFERENCE = [
    # (url, destination, gunzip)
    ("https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz",
     "data/external/s288c/sacCer3.fa", True),
    ("http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/saccharomyces_cerevisiae.gff.gz",
     "data/external/sgd/saccharomyces_cerevisiae.gff.gz", False),
    ("https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0030-5/MediaObjects/"
     "41586_2018_30_MOESM3_ESM.xls", "data/41586_2018_30_MOESM3_ESM.xls", False),
    ("https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-025-09637-0/MediaObjects/"
     "41586_2025_9637_MOESM3_ESM.xlsx", "data/41586_2025_9637_MOESM3_ESM.xlsx", False),
    # Helsen et al. 2026 Supplementary Data 4 (CEN annotations, used by B and F2) and 6 (assemblies)
    ("https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-025-09779-1/MediaObjects/"
     "41586_2025_9779_MOESM5_ESM.xlsx", "data/external/helsen2026/SD4_scer_centromeres.xlsx", False),
    ("https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-025-09779-1/MediaObjects/"
     "41586_2025_9779_MOESM7_ESM.xlsx", "data/external/helsen2026/SD6_scer_assemblies.xlsx", False),
]
LOEGLER_TABLE = "data/41586_2025_9637_MOESM3_ESM.xlsx"

SCRAP_RE = re.compile(r"(?:^|/)([^/]+)\.asm01\.(collapsed|HP0|HP1|HP2)\.nuclear_genome\.tidy\.fa\.gz$")
CHUNK = 1 << 20


# ---------------------------------------------------------------- selection
def scrap_targets(names: list[str], haplotypes: str = "one") -> dict[str, dict]:
    """Archive member -> {sample, assembly, name}. One per strain: collapsed > HP0."""
    by_strain: dict[str, dict[str, str]] = {}
    for m in names:
        hit = SCRAP_RE.search(m)
        if hit:
            by_strain.setdefault(hit.group(1), {})[hit.group(2)] = m
    out = {}
    for name, haps in sorted(by_strain.items()):
        keep = []
        for h in ("collapsed", "HP0"):
            if h in haps:
                keep.append(h)
                break
        if haplotypes == "all":
            keep += [h for h in ("HP1", "HP2") if h in haps]
        for h in keep:
            sample = name if h == "HP0" else f"{name}_{h}"
            out[haps[h]] = {"sample": sample, "assembly": sample, "name": name, "primary": h == keep[0]}
    return out


def loegler_targets(table: pd.DataFrame, haplotypes: str = "one") -> dict[str, dict]:
    """Final-assembly file name -> {sample, assembly, name}, from Loegler Table S2."""
    out = {}
    for _, r in table.iterrows():
        name, hap = str(r["Standardized_name"]).strip(), r["Haplotype"]
        fn = str(r["Filename"]).strip()
        if pd.isna(hap):
            out[fn] = {"sample": name, "assembly": name, "name": name, "primary": True}
        elif hap == "HP1" or haplotypes == "all":
            sample = f"{name}_{hap}"
            out[fn] = {"sample": sample, "assembly": sample, "name": name, "primary": hap == "HP1"}
    return out


def read_loegler_table(path: Path) -> pd.DataFrame:
    t = pd.read_excel(path, sheet_name="TableS2", header=2)
    return t[["Standardized_name", "Haplotype", "Filename"]].dropna(subset=["Filename"])


# ScRAP names whose `_1b`/`_3`-style suffix is needed to keep strains unique (checked on the
# 142-strain member list, 2026-11-16)
SCRAP_KEEP_SUFFIX = {"ALH", "CGH"}


def strain_names(targets: dict[str, dict], dataset: str) -> dict[str, str]:
    """sample -> strain (a unique key). ScRAP: name without a `_1a`-style suffix (`CCC_1a` ->
    `CCC`, the Peter 2018 key) except SCRAP_KEEP_SUFFIX. Loegler: Standardized_name as is (it has
    both `AGA` and `AGA_1a`). Extra haplotypes (`--haplotypes all`) keep their sample name."""
    def short(n):
        b = re.sub(r"_\d+[a-z]*$", "", n)
        return b if dataset == "scrap142" and b not in SCRAP_KEEP_SUFFIX else n
    return {t["sample"]: (short(t["name"]) if t["primary"] else t["sample"]) for t in targets.values()}


def apply_filter(targets: dict[str, dict], strains: list[str] | None) -> dict[str, dict]:
    if not strains:
        return targets
    want = set(strains)
    out = {m: t for m, t in targets.items() if t["name"] in want or t["sample"] in want}
    missing = want - {t["name"] for t in out.values()} - {t["sample"] for t in out.values()}
    if missing:
        sys.exit(f"Not in this set: {', '.join(sorted(missing))}")
    return out


# ---------------------------------------------------------------- streaming
class HashingReader:
    """File-like wrapper that md5-hashes everything read through it."""

    def __init__(self, raw):
        self.raw, self.md5, self.n = raw, hashlib.md5(), 0

    def read(self, size=-1):
        b = self.raw.read(size)
        self.md5.update(b)
        self.n += len(b)
        return b


def zenodo_file(record: str, key: str) -> tuple[str, str, int]:
    with urllib.request.urlopen(f"https://zenodo.org/api/records/{record}", timeout=60) as r:
        meta = json.load(r)
    f = next(x for x in meta["files"] if x["key"] == key)
    return f["links"]["self"], f["checksum"].removeprefix("md5:"), int(f["size"])


def write_member(src, dest: Path, already_gz: bool) -> None:
    tmp = dest.with_name(dest.name + ".tmp")
    with open(tmp, "wb") as out:
        if already_gz:
            shutil.copyfileobj(src, out, CHUNK)
        else:
            with gzip.GzipFile(fileobj=out, mode="wb", compresslevel=6, mtime=0) as gz:
                shutil.copyfileobj(src, gz, CHUNK)
    tmp.replace(dest)


def extract(stream, targets: dict[str, dict], key, outdir: Path, limit: int | None) -> tuple[list[dict], bool]:
    """Stream a tar.gz, write wanted members; returns (written-or-present targets, read to the end)."""
    outdir.mkdir(parents=True, exist_ok=True)
    done, todo = [], dict(targets)
    for t in list(todo.values()):
        if (outdir / f"{t['sample']}.fasta.gz").exists():
            done.append(t)
    todo = {m: t for m, t in todo.items() if t not in done}
    if limit is not None:
        if len(done) >= limit:
            return done[:limit], False
    if not todo:
        return done, False
    with tarfile.open(fileobj=stream, mode="r|gz") as tar:
        for member in tar:
            k = key(member.name)
            if not member.isfile() or k not in todo:
                continue
            t = todo.pop(k)
            dest = outdir / f"{t['sample']}.fasta.gz"
            write_member(tar.extractfile(member), dest, member.name.endswith(".gz"))
            done.append(t)
            print(f"  {len(done):5d}  {t['sample']}  ({member.size / 1e6:.1f} MB)", flush=True)
            if not todo or (limit is not None and len(done) >= limit):
                return done, False
    if todo:
        print(f"WARNING: {len(todo)} wanted assemblies not found in the archive: "
              f"{', '.join(sorted(t['sample'] for t in todo.values())[:10])} ...", file=sys.stderr)
    return done, True


def download_assemblies(args) -> None:
    record, fname = ZENODO[args.set]
    if args.set == "scrap142":
        listing = None
        if args.archive:
            with tarfile.open(args.archive, "r:gz") as tar:
                listing = tar.getnames()
        targets = scrap_targets(listing, args.haplotypes) if listing else None
        key = lambda name: name  # noqa: E731
    else:
        tab = REPO / LOEGLER_TABLE
        if not tab.exists():
            sys.exit(f"{tab} missing: run `--set reference` first")
        targets = loegler_targets(read_loegler_table(tab), args.haplotypes)
        key = lambda name: Path(name).name if "/Assemblies_SecondBest" not in name else None  # noqa: E731

    outdir = Path(args.out or REPO / f"data/assemblies_{args.set}")
    if targets is None:
        # ScRAP member list is only known from the archive: stream once, select on the fly
        targets = _scrap_stream(args, record, fname, outdir)
    else:
        targets = apply_filter(targets, args.strains)
        print(f"{args.set}: {len(targets)} assemblies wanted -> {outdir}")
        if args.archive:
            with open(args.archive, "rb") as fh:
                done, full = extract(fh, targets, key, outdir, args.limit)
        else:
            url, md5, size = zenodo_file(record, fname)
            print(f"streaming {fname} ({size / 1e9:.1f} GB) from Zenodo {record}")
            with urllib.request.urlopen(url, timeout=120) as r:
                hr = HashingReader(r)
                done, full = extract(hr, targets, key, outdir, args.limit)
                if full:
                    shutil.copyfileobj(hr, open("/dev/null", "wb"), CHUNK)  # tar end padding
                    ok = hr.md5.hexdigest() == md5
                    print(f"md5 {'OK' if ok else 'MISMATCH'} ({hr.md5.hexdigest()})")
                    if not ok:
                        sys.exit("archive md5 mismatch: delete the outputs and download again")
        targets = {t["sample"]: t for t in done}
    write_samples(targets, args.set)


def _scrap_stream(args, record, fname, outdir) -> dict:
    """ScRAP without a local archive: decide per member while streaming (members are named
    <strain>.asm01.<hap>...)."""
    url, md5, size = zenodo_file(record, fname)
    print(f"streaming {fname} ({size / 1e9:.1f} GB) from Zenodo {record} -> {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)
    want = set(args.strains) if args.strains else None
    seen: dict[str, dict[str, str]] = {}
    with urllib.request.urlopen(url, timeout=120) as r:
        hr = HashingReader(r)
        with tarfile.open(fileobj=hr, mode="r|gz") as tar:
            for member in tar:
                hit = SCRAP_RE.search(member.name)
                if not member.isfile() or not hit:
                    continue
                name, hap = hit.groups()
                if want and name not in want:
                    continue
                seen.setdefault(name, {})[hap] = member.name
                if hap not in ("collapsed", "HP0") and args.haplotypes != "all":
                    continue
                sample = name if hap == "HP0" else f"{name}_{hap}"
                dest = outdir / f"{sample}.fasta.gz"
                if not dest.exists():
                    write_member(tar.extractfile(member), dest, True)
                    print(f"  {sample}", flush=True)
                n = len(scrap_targets([m for h in seen.values() for m in h.values()], args.haplotypes))
                if args.limit is not None and n >= args.limit:
                    break
            else:
                shutil.copyfileobj(hr, open("/dev/null", "wb"), CHUNK)
                ok = hr.md5.hexdigest() == md5
                print(f"md5 {'OK' if ok else 'MISMATCH'} ({hr.md5.hexdigest()})")
                if not ok:
                    sys.exit("archive md5 mismatch: delete the outputs and download again")
    # ScRAP strains have either `collapsed` (+HP/HP1/HP2) or `HP0`, never both
    targets = scrap_targets([m for h in seen.values() for m in h.values()], args.haplotypes)
    if args.limit is not None:
        targets = dict(list(targets.items())[:args.limit])
    return {t["sample"]: t for t in targets.values()}


def write_samples(targets: dict[str, dict], name: str) -> None:
    strains = strain_names(targets, name)
    rows = [{"accession": t["sample"], "assembly": t["assembly"], "strain": strains[t["sample"]]}
            for t in sorted(targets.values(), key=lambda t: t["sample"])]
    path = REPO / f"data/samples_{name}.tsv"
    pd.DataFrame(rows, columns=["accession", "assembly", "strain"]).to_csv(path, sep="\t", index=False)
    print(f"{len(rows)} samples -> {path.relative_to(REPO)}")


# ---------------------------------------------------------------- reference
def download_reference() -> None:
    for url, dest, gunzip in REFERENCE:
        dest = REPO / dest
        if dest.exists():
            print(f"skip {dest.relative_to(REPO)} (exists)")
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        tmp = dest.with_name(dest.name + ".tmp")
        print(f"{url} -> {dest.relative_to(REPO)}")
        with urllib.request.urlopen(url, timeout=120) as r, open(tmp, "wb") as out:
            shutil.copyfileobj(gzip.GzipFile(fileobj=r) if gunzip else r, out, CHUNK)
        tmp.replace(dest)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--set", required=True, choices=["reference", *ZENODO])
    ap.add_argument("--out", help="assembly directory (default data/assemblies_<set>)")
    ap.add_argument("--archive", help="local copy of the Zenodo tar.gz instead of streaming")
    ap.add_argument("--limit", type=int, help="stop after N assemblies (test runs)")
    ap.add_argument("--strains", type=lambda s: [x for x in s.split(",") if x],
                    help="comma-separated strain names (e.g. ADS,CCC_1a)")
    ap.add_argument("--haplotypes", choices=["one", "all"], default="one")
    args = ap.parse_args()
    if args.set == "reference":
        download_reference()
    else:
        download_assemblies(args)


if __name__ == "__main__":
    main()
