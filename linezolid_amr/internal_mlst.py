"""In-house MLST implementation backed by PubMLST schemes.

Replaces the optional external `mlst` (Seemann) dependency to avoid the
samtools 0.1.x dependency conflict that prevents distributing them together
in the same conda env.

Architecture
============

For each supported organism we ship a small **scheme directory** under
``linezolid_amr/data/mlst_schemes/<scheme_name>/`` containing:

  - ``profiles.tsv``        — ST profile table (ST<TAB>locus1<TAB>locus2<TAB>...)
  - ``loci.txt``            — one locus name per line, defines column order
  - ``<locus>.fasta.gz``    — gzipped FASTA of all alleles for that locus,
                              headers like ``>arcC_1``, ``>arcC_2``, ...

The runtime calls ``blastn`` to align the assembly against the per-locus FASTA,
picks the best hit per locus, and looks the resulting allele tuple up in
``profiles.tsv``. The output is a :class:`MlstResult` dataclass that mirrors
what the legacy ``mlst.py`` wrapper produced.

PubMLST is the canonical source. The schemes are populated by
``fetch_pubmlst_schemes()`` (also exposed as a CLI subcommand) and bundled
inside the package so installs work entirely offline.
"""

from __future__ import annotations

import csv
import gzip
import io
import json
import os
import shutil
import subprocess
import tempfile
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass
from importlib import resources
from pathlib import Path
from typing import Iterable

# Mapping: AMRFinderPlus-style organism → PubMLST seqdef database + short scheme name
PUBMLST_DBS = {
    "Staphylococcus_aureus":    {"db": "pubmlst_saureus_seqdef",    "scheme": 1, "name": "saureus"},
    "Enterococcus_faecalis":    {"db": "pubmlst_efaecalis_seqdef",  "scheme": 1, "name": "efaecalis"},
    "Enterococcus_faecium":     {"db": "pubmlst_efaecium_seqdef",   "scheme": 1, "name": "efaecium"},
    "Streptococcus_pneumoniae": {"db": "pubmlst_spneumoniae_seqdef", "scheme": 1, "name": "spneumoniae"},
}

PUBMLST_BASE = "https://rest.pubmlst.org"
USER_AGENT = "linezolid-amr/0.1 (github.com/iowa69/linezolid-amr)"

PACKAGE_SCHEMES_RESOURCE = "linezolid_amr.data.mlst_schemes"

# Match thresholds — aligned 1:1 with tseemann/mlst defaults so calls are
# directly comparable. See https://github.com/tseemann/mlst/blob/master/bin/mlst
MIN_IDENTITY = 95.0   # --minid 95 (% DNA identity threshold)
MIN_COVERAGE = 50.0   # --mincov 50 (% coverage threshold)
EXACT_IDENTITY = 100.0
EXACT_COVERAGE = 100.0


# ----------------------------------------------------------------------
# Public dataclass — keeps the same shape as legacy mlst.MlstResult
# ----------------------------------------------------------------------

@dataclass
class MlstResult:
    file: str
    scheme: str
    st: str
    alleles: dict[str, str]
    raw: str

    @property
    def organism(self) -> str | None:
        """AMRFinderPlus organism name corresponding to this scheme."""
        for org, info in PUBMLST_DBS.items():
            if info["name"] == self.scheme:
                return org
        return None


class MlstNoMatch(ValueError):
    """No allele matched well enough at any locus."""


class MlstUnsupportedOrganism(ValueError):
    """User asked for MLST on an organism outside the LZD-supported list."""


# ----------------------------------------------------------------------
# Scheme location
# ----------------------------------------------------------------------

def schemes_override_dir() -> Path | None:
    """Return override schemes dir if user set ``LINEZOLID_AMR_MLST_DIR``."""
    env = os.environ.get("LINEZOLID_AMR_MLST_DIR")
    if env:
        return Path(env).expanduser().resolve()
    return None


def scheme_dir(scheme_name: str) -> Path:
    """Resolve the directory containing scheme files for ``scheme_name``.

    Lookup order:
      1. ``$LINEZOLID_AMR_MLST_DIR/<scheme>/``  (user override)
      2. bundled package resource

    Returns a filesystem Path. Does not check existence.
    """
    override = schemes_override_dir()
    if override is not None:
        candidate = override / scheme_name
        if candidate.exists():
            return candidate
    return Path(str(resources.files(PACKAGE_SCHEMES_RESOURCE).joinpath(scheme_name)))


def scheme_available(scheme_name: str) -> bool:
    d = scheme_dir(scheme_name)
    return d.exists() and (d / "profiles.tsv").exists() and (d / "loci.txt").exists()


def all_schemes_available() -> bool:
    return all(scheme_available(info["name"]) for info in PUBMLST_DBS.values())


# ----------------------------------------------------------------------
# PubMLST fetcher
# ----------------------------------------------------------------------

def _http_get_bytes(url: str, retries: int = 3, backoff: float = 1.5) -> bytes:
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=60) as resp:
                return resp.read()
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(backoff * (attempt + 1))
    raise RuntimeError(f"PubMLST request failed after {retries} retries: {url}\n  {last}")


def fetch_pubmlst_scheme(organism: str, out_root: Path, force: bool = False) -> dict:
    """Fetch one MLST scheme from PubMLST. Returns a summary dict."""
    if organism not in PUBMLST_DBS:
        raise MlstUnsupportedOrganism(
            f"Organism '{organism}' is not in the supported list "
            f"({', '.join(PUBMLST_DBS)})."
        )
    info = PUBMLST_DBS[organism]
    scheme_name = info["name"]
    out = out_root / scheme_name
    out.mkdir(parents=True, exist_ok=True)

    if not force and (out / "profiles.tsv").exists() and (out / "loci.txt").exists():
        return {"organism": organism, "status": "cached", "dir": str(out)}

    # 1) Scheme JSON to get locus list
    sj = json.loads(_http_get_bytes(f"{PUBMLST_BASE}/db/{info['db']}/schemes/{info['scheme']}").decode())
    locus_urls = sj["loci"]
    locus_names = [u.rsplit("/", 1)[-1] for u in locus_urls]
    (out / "loci.txt").write_text("\n".join(locus_names) + "\n")

    # 2) Profile CSV → TSV
    csv_bytes = _http_get_bytes(sj["profiles_csv"])
    text = csv_bytes.decode("utf-8")
    # PubMLST returns tab-separated despite the "csv" name — defensively split on any whitespace tab
    # We re-emit as a clean TSV with only ST + the loci columns.
    lines = [ln for ln in text.splitlines() if ln.strip()]
    header = lines[0].split("\t")
    # ST column is always present, sometimes capitalised. Find it.
    if "ST" in header:
        st_idx = header.index("ST")
    else:
        st_idx = 0  # fallback: first column
    locus_idx = [header.index(l) for l in locus_names]
    with (out / "profiles.tsv").open("w") as fh:
        fh.write("ST\t" + "\t".join(locus_names) + "\n")
        for ln in lines[1:]:
            parts = ln.split("\t")
            try:
                row = [parts[st_idx]] + [parts[i] for i in locus_idx]
            except IndexError:
                continue  # malformed line
            fh.write("\t".join(row) + "\n")

    # 3) Per-locus allele FASTAs (gzipped)
    n_alleles = 0
    for locus, locus_url in zip(locus_names, locus_urls):
        fa = _http_get_bytes(f"{locus_url}/alleles_fasta")
        with gzip.open(out / f"{locus}.fasta.gz", "wb") as fh:
            fh.write(fa)
        # crude allele count
        n_alleles += fa.count(b">")
        time.sleep(0.2)  # be gentle to PubMLST

    return {
        "organism": organism,
        "status": "fetched",
        "dir": str(out),
        "loci": locus_names,
        "n_alleles": n_alleles,
    }


def fetch_pubmlst_schemes(
    organisms: Iterable[str] | None = None,
    out_root: Path | None = None,
    force: bool = False,
) -> list[dict]:
    """Fetch schemes for all (or a subset) of supported organisms.

    If ``out_root`` is None, writes into the user override dir (creating it).
    """
    if out_root is None:
        env = os.environ.get("LINEZOLID_AMR_MLST_DIR")
        out_root = Path(env).expanduser() if env else Path.home() / ".local" / "share" / "linezolid-amr" / "mlst_schemes"
    out_root.mkdir(parents=True, exist_ok=True)

    targets = list(organisms) if organisms else list(PUBMLST_DBS)
    results: list[dict] = []
    for org in targets:
        try:
            results.append(fetch_pubmlst_scheme(org, out_root, force=force))
        except Exception as e:  # noqa: BLE001
            results.append({"organism": org, "status": "error", "error": str(e)})
    return results


# ----------------------------------------------------------------------
# BLAST runner & MLST caller
# ----------------------------------------------------------------------

def blastn_available() -> bool:
    return shutil.which("blastn") is not None and shutil.which("makeblastdb") is not None


def _ensure_blast_db(allele_fasta_gz: Path, work: Path) -> Path:
    """Decompress an allele FASTA and build a BLAST nucleotide database."""
    base = allele_fasta_gz.stem  # strip .gz → .fasta
    db_prefix = work / Path(base).stem  # strip .fasta
    fasta = work / base
    if not fasta.exists():
        with gzip.open(allele_fasta_gz, "rb") as gin, fasta.open("wb") as fout:
            shutil.copyfileobj(gin, fout)
    # makeblastdb is idempotent — silently skip if .nhr already exists
    if not (db_prefix.with_suffix(".nhr").exists()):
        subprocess.run(
            ["makeblastdb", "-in", str(fasta), "-dbtype", "nucl", "-out", str(db_prefix)],
            check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
        )
    return db_prefix


def _best_allele(assembly: Path, blast_db: Path, locus: str, threads: int) -> list[tuple[str, float, float]]:
    """Run blastn(assembly vs allele DB). Return all qualifying hits sorted best-first.

    Each hit is ``(allele_id, identity, coverage_pct)``.
    BLAST options match tseemann/mlst defaults verbatim:
      -ungapped -dust no -word_size 32 -evalue 1E-20 -perc_identity 95
    """
    cmd = [
        "blastn",
        "-query", str(assembly),
        "-db", str(blast_db),
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen",
        "-perc_identity", str(MIN_IDENTITY),
        "-num_threads", str(threads),
        "-ungapped",
        "-dust", "no",
        "-word_size", "32",
        "-evalue", "1E-20",
        "-max_target_seqs", "100000",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=True)
    hits: list[tuple[str, float, float]] = []
    for line in proc.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 13:
            continue
        sseqid = parts[1]
        pident = float(parts[2])
        aln_len = int(parts[3])
        slen = int(parts[12])
        coverage = aln_len / slen * 100 if slen else 0.0
        if coverage < MIN_COVERAGE:
            continue
        allele_id = sseqid.split("_")[-1]
        hits.append((allele_id, pident, coverage))
    # Sort best-first by (coverage desc, identity desc, allele number asc)
    hits.sort(key=lambda h: (-h[2], -h[1], int(h[0]) if h[0].isdigit() else 10**9))
    return hits


def _annotate_allele(hits: list[tuple[str, float, float]]) -> str:
    """Map BLAST hits to tseemann/mlst-style locus notation.

    - ``n``         : a single perfect (100% id + 100% cov) match
    - ``n,m``       : two or more perfect matches at different allele numbers
    - ``n?``        : single full-identity but partial-coverage hit
    - ``~n``        : best near-perfect (≥95% id, ≥50% cov) but not perfect
    - ``-``         : no qualifying hit
    """
    if not hits:
        return "-"
    perfects = [h for h in hits if h[1] >= EXACT_IDENTITY and h[2] >= EXACT_COVERAGE]
    if perfects:
        uniq = sorted({h[0] for h in perfects}, key=lambda x: int(x) if x.isdigit() else 10**9)
        return ",".join(uniq) if len(uniq) > 1 else uniq[0]
    # full identity but partial coverage → "n?"
    partial_cov = [h for h in hits if h[1] >= EXACT_IDENTITY]
    if partial_cov:
        return f"{partial_cov[0][0]}?"
    # near-perfect (identity 95-100 OR coverage 50-100) → "~n"
    return f"~{hits[0][0]}"


def _load_profiles(scheme_path: Path) -> tuple[list[str], dict[tuple[str, ...], str]]:
    """Return (locus_order, {allele_tuple -> ST})."""
    loci = (scheme_path / "loci.txt").read_text().splitlines()
    loci = [l.strip() for l in loci if l.strip()]
    profiles: dict[tuple[str, ...], str] = {}
    with (scheme_path / "profiles.tsv").open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            try:
                tup = tuple(str(row[l]) for l in loci)
            except KeyError:
                continue
            profiles[tup] = str(row["ST"])
    return loci, profiles


def run_internal_mlst(assembly: Path, organism: str, threads: int = 4) -> MlstResult:
    """Run our PubMLST-based MLST. Returns an MlstResult."""
    # Validate organism FIRST so callers get the right exception type, even
    # if blastn happens to be missing on the system.
    if organism not in PUBMLST_DBS:
        raise MlstUnsupportedOrganism(
            f"Organism '{organism}' is outside the supported MLST set "
            f"({', '.join(PUBMLST_DBS)}). Pass --organism manually if you only need AMRFinderPlus."
        )
    if not blastn_available():
        raise RuntimeError("blastn / makeblastdb not found on PATH. Install via `conda install -c bioconda blast`.")
    info = PUBMLST_DBS[organism]
    scheme_name = info["name"]
    sdir = scheme_dir(scheme_name)
    if not scheme_available(scheme_name):
        raise FileNotFoundError(
            f"MLST scheme '{scheme_name}' not bundled in this build and not present in "
            f"$LINEZOLID_AMR_MLST_DIR. Run `linezolid-amr fetch-mlst-schemes` first."
        )

    locus_order, profiles = _load_profiles(sdir)

    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        alleles_found: dict[str, str] = {}
        for locus in locus_order:
            allele_fasta_gz = sdir / f"{locus}.fasta.gz"
            if not allele_fasta_gz.exists():
                alleles_found[locus] = "-"
                continue
            blast_db = _ensure_blast_db(allele_fasta_gz, work)
            hits = _best_allele(assembly, blast_db, locus, threads)
            alleles_found[locus] = _annotate_allele(hits)

    tup = tuple(alleles_found[l] for l in locus_order)
    st, st_note = _assign_st(tup, locus_order, profiles)

    raw = "\t".join([str(assembly), scheme_name, st] + [f"{l}({alleles_found[l]})" for l in locus_order])
    return MlstResult(
        file=str(assembly),
        scheme=scheme_name,
        st=st,
        alleles=alleles_found,
        raw=raw,
    )


def _assign_st(
    tup: tuple[str, ...],
    locus_order: list[str],
    profiles: dict[tuple[str, ...], str],
) -> tuple[str, str]:
    """Decide the ST given the per-locus call tuple.

    Exact integer-tuple lookup in the PubMLST profile table, with one extra
    rule that ``mlst`` (Seemann) does NOT implement: PubMLST encodes a
    genuinely-deleted locus as allele ``0`` in the profile table (e.g.
    *E. faecium* ST1478 has ``pstS=0``). Our locus caller emits ``-`` when no
    BLAST hit clears threshold — which is the read-level expectation of a
    deleted gene. We therefore retry the lookup with ``-`` substituted by
    ``0`` so null-allele STs like 1478/1421/1422/... get assigned correctly.

    Partial matches (~n), ambiguous calls (n?) still yield ``-``.
    """
    # First pass: strict exact match (no missing loci).
    if all(a.isdigit() for a in tup) and tup in profiles:
        return profiles[tup], ""
    # Second pass: allow '-' (missing locus) to match PubMLST null-allele '0'.
    # Only kicks in if every non-missing locus is a clean integer; otherwise we
    # would risk mapping noise into a real ST.
    if all(a.isdigit() or a == "-" for a in tup) and any(a == "-" for a in tup):
        normalized = tuple("0" if a == "-" else a for a in tup)
        if normalized in profiles:
            return profiles[normalized], ""
    return "-", ""


def infer_organism(assembly: Path, threads: int = 4) -> tuple[str, MlstResult]:
    """Try every supported organism's MLST scheme, pick the one with the most
    exact allele hits. If none have any hit, raise MlstNoMatch."""
    best: tuple[str, MlstResult] | None = None
    best_score = -1
    for organism in PUBMLST_DBS:
        try:
            res = run_internal_mlst(assembly, organism, threads=threads)
        except Exception:
            continue
        # score: number of perfect-allele matches (digit-only IDs)
        score = sum(1 for v in res.alleles.values() if v.isdigit())
        # tie-break: a real ST beats "-"
        if score > best_score or (score == best_score and best and res.st != "-" and best[1].st == "-"):
            best = (organism, res)
            best_score = score
    if not best or best_score == 0:
        raise MlstNoMatch(
            "No MLST scheme matched the assembly above thresholds. "
            "Either the species is outside our supported list or the assembly is fragmented. "
            "Pass --organism explicitly to bypass inference."
        )
    return best
