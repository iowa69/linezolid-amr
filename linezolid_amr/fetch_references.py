"""Fetch and prepare per-species 23S rRNA reference sequences from NCBI.

For each supported organism we:

  1) efetch the verified 23S rRNA slice from its reference genome (NCBI E-utils)
  2) efetch the E. coli K-12 rrlB 23S to act as the coordinate master
  3) globally align species-23S vs E. coli-23S to derive a position map
  4) write:
       <cache>/<organism>_23S.fasta
       <cache>/<organism>_23S_position_map.tsv  (ecoli_pos -> species_pos)
       <cache>/<organism>_23S_lzd_positions.bed (species coords of LZD-resistance sites)

The position mapping is what lets us report mutations in the clinically standard
E. coli 23S numbering even though reads are aligned to species-specific references.
"""

from __future__ import annotations

import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Iterable

from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from linezolid_amr.references import (
    EcoliRef,
    OrganismRef,
    cache_dir,
    ecoli_fasta_path,
    get_ecoli_reference,
    get_linezolid_positions,
    get_organism,
    list_organisms,
    organism_bed_path,
    organism_fasta_path,
    organism_position_map_path,
)

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
USER_AGENT = "linezolid-amr/0.1 (github.com/iowa69/linezolid-amr)"


def _http_get(url: str, retries: int = 3, backoff: float = 2.0) -> bytes:
    last_err = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=60) as resp:
                return resp.read()
        except Exception as e:  # noqa: BLE001
            last_err = e
            time.sleep(backoff * (attempt + 1))
    raise RuntimeError(f"NCBI request failed after {retries} retries: {url}\n  {last_err}")


def efetch_fasta_slice(
    accession: str,
    start: int,
    end: int,
    strand: str,
    api_key: str | None = None,
) -> SeqRecord:
    """Fetch a FASTA slice from NCBI nuccore."""
    params = {
        "db": "nuccore",
        "id": accession,
        "rettype": "fasta",
        "retmode": "text",
        "seq_start": str(start),
        "seq_stop": str(end),
        "strand": "1" if strand == "+" else "2",
    }
    if api_key:
        params["api_key"] = api_key
    import io
    url = f"{EUTILS_BASE}/efetch.fcgi?{urllib.parse.urlencode(params)}"
    data = _http_get(url)
    text = data.decode("utf-8")
    records = list(SeqIO.parse(io.StringIO(text), "fasta"))
    if not records:
        raise RuntimeError(f"No FASTA returned for {accession}:{start}-{end} strand {strand}")
    return records[0]


def _pairwise_align(query: Seq, target: Seq) -> tuple[list[int], list[int]]:
    """Global alignment of query vs target. Returns (q_to_t, t_to_q) 0-based position maps.

    For each q index i, q_to_t[i] is the 0-based target index aligned to it (or -1 for a gap).
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -1
    aln = aligner.align(query, target)[0]
    # Walk alignment coordinates
    q_to_t = [-1] * len(query)
    t_to_q = [-1] * len(target)
    qi, ti = 0, 0
    # aln.aligned -> ((q_blocks), (t_blocks)) of equal-length tuple-of-(start,end)
    q_blocks, t_blocks = aln.aligned
    for (qs, qe), (ts, te) in zip(q_blocks, t_blocks):
        block = qe - qs
        for k in range(block):
            q_to_t[qs + k] = ts + k
            t_to_q[ts + k] = qs + k
    return q_to_t, t_to_q


def _write_fasta(record: SeqRecord, path: Path, label: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    record = SeqRecord(record.seq, id=label, description="")
    with path.open("w") as fh:
        SeqIO.write([record], fh, "fasta")


def fetch_ecoli_reference(api_key: str | None = None, force: bool = False) -> Path:
    ref = get_ecoli_reference()
    out = ecoli_fasta_path()
    if out.exists() and not force:
        return out
    rec = efetch_fasta_slice(ref.genome_accession, ref.start, ref.end, ref.strand, api_key=api_key)
    if abs(len(rec.seq) - ref.expected_length_bp) > 20:
        print(
            f"warning: E. coli 23S length {len(rec.seq)} differs from expected "
            f"{ref.expected_length_bp}",
            file=sys.stderr,
        )
    _write_fasta(rec, out, ref.label)
    return out


def fetch_organism_reference(
    organism: str, ecoli_seq: Seq, api_key: str | None = None, force: bool = False
) -> dict:
    org: OrganismRef = get_organism(organism)
    fasta_path = organism_fasta_path(organism)
    bed_path = organism_bed_path(organism)
    map_path = organism_position_map_path(organism)
    if fasta_path.exists() and bed_path.exists() and map_path.exists() and not force:
        return {"organism": organism, "status": "cached", "fasta": str(fasta_path)}

    rec = efetch_fasta_slice(
        org.genome_accession, org.start, org.end, org.strand, api_key=api_key
    )
    if abs(len(rec.seq) - org.expected_length_bp) > 30:
        print(
            f"warning: {organism} 23S length {len(rec.seq)} differs from expected "
            f"{org.expected_length_bp}",
            file=sys.stderr,
        )

    label = f"{organism}_23S"
    _write_fasta(rec, fasta_path, label)

    # Align species 23S (query) vs E. coli 23S (target) to derive position map
    species_seq = rec.seq.upper()
    target = ecoli_seq.upper()
    q_to_t, t_to_q = _pairwise_align(species_seq, target)

    # E. coli numbering is 1-based; t_to_q is 0-based target->query.
    map_path.parent.mkdir(parents=True, exist_ok=True)
    with map_path.open("w") as fh:
        fh.write("ecoli_position\tspecies_position\tspecies_base\tecoli_base\n")
        for ti in range(len(target)):
            qi = t_to_q[ti]
            if qi < 0:
                fh.write(f"{ti + 1}\t.\t-\t{target[ti]}\n")
            else:
                fh.write(f"{ti + 1}\t{qi + 1}\t{species_seq[qi]}\t{target[ti]}\n")

    # Build BED of canonical LZD positions in species coordinates
    positions = get_linezolid_positions()
    with bed_path.open("w") as fh:
        fh.write("#chrom\tstart\tend\tname\tref_base\tresistance_bases\tdrug\n")
        for p in positions:
            ti = p.ecoli_position - 1  # 0-based
            if ti < 0 or ti >= len(target):
                continue
            qi = t_to_q[ti]
            if qi < 0:
                continue
            species_pos_1b = qi + 1
            name = f"23S_E{p.ecoli_position}_{p.ref_base}_to_{'/'.join(p.resistance_bases)}"
            fh.write(
                f"{label}\t{qi}\t{species_pos_1b}\t{name}\t{p.ref_base}\t"
                f"{','.join(p.resistance_bases)}\t{p.drug}\n"
            )

    return {
        "organism": organism,
        "status": "fetched",
        "fasta": str(fasta_path),
        "bed": str(bed_path),
        "position_map": str(map_path),
        "length_bp": len(species_seq),
    }


def fetch_all(
    organisms: Iterable[str] | None = None,
    api_key: str | None = None,
    force: bool = False,
) -> list[dict]:
    """Fetch all (or a subset of) supported organisms. Returns list of result dicts."""
    targets = list(organisms) if organisms else list_organisms()
    out: list[dict] = []

    ecoli_fa = fetch_ecoli_reference(api_key=api_key, force=force)
    ecoli_rec = next(SeqIO.parse(str(ecoli_fa), "fasta"))
    ecoli_seq = ecoli_rec.seq

    for organism in targets:
        try:
            res = fetch_organism_reference(organism, ecoli_seq, api_key=api_key, force=force)
        except Exception as e:  # noqa: BLE001
            res = {"organism": organism, "status": "error", "error": str(e)}
        out.append(res)
        # be polite to NCBI between accessions
        time.sleep(0.4)
    return out
