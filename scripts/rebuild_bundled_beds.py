#!/usr/bin/env python3
"""Regenerate the bundled per-species LZD target BEDs from loci.json.

Runs entirely offline: the species coordinate for each E. coli position is read
from the bundled position maps (which `fetch-references` produced by global
alignment), and the reference base is read from the species' own 23S FASTA.

Use this after editing the curated position list in loci.json. Network access is
only needed if the underlying references themselves change, which is what
`linezolid-amr fetch-references` is for.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from linezolid_amr import references as ref_mod  # noqa: E402

HEADER = "#chrom\tstart\tend\tname\tref_base\tresistance_bases\tdrug\tevidence_tier\n"


def read_fasta(path: Path) -> tuple[str, str]:
    lines = path.read_text().splitlines()
    name = lines[0][1:].split()[0]
    return name, "".join(l.strip() for l in lines[1:] if not l.startswith(">")).upper()


def load_map(path: Path) -> dict[int, int]:
    """E. coli position -> species position (1-based); skips gaps."""
    out: dict[int, int] = {}
    for line in path.read_text().splitlines()[1:]:
        e_pos, sp_pos, _sp_base, _e_base = line.split("\t")
        if sp_pos not in (".", ""):
            out[int(e_pos)] = int(sp_pos)
    return out


def main() -> int:
    positions = ref_mod.get_linezolid_positions()
    ecoli_seq = read_fasta(ref_mod.ecoli_fasta_path())[1]

    # Guard: the curated E. coli bases must match the E. coli reference, or the
    # whole coordinate framework is untrustworthy.
    bad = [
        f"{p.ecoli_position}: curated {p.ref_base}, reference {ecoli_seq[p.ecoli_position - 1]}"
        for p in positions
        if ecoli_seq[p.ecoli_position - 1] != p.ref_base.upper()
    ]
    if bad:
        print("ABORT — loci.json disagrees with the E. coli reference:", file=sys.stderr)
        for b in bad:
            print("  " + b, file=sys.stderr)
        return 1

    for organism in ref_mod.list_organisms():
        contig, species_seq = read_fasta(ref_mod.organism_fasta_path(organism))
        emap = load_map(ref_mod.organism_position_map_path(organism))
        bed_path = ref_mod.organism_bed_path(organism)

        rows, skipped = [], []
        for p in sorted(positions, key=lambda x: x.ecoli_position):
            sp = emap.get(p.ecoli_position)
            if sp is None:
                skipped.append(p.ecoli_position)
                continue
            species_base = species_seq[sp - 1]
            name = f"23S_E{p.ecoli_position}_{p.ref_base}_to_{'/'.join(p.resistance_bases)}"
            rows.append(
                f"{contig}\t{sp - 1}\t{sp}\t{name}\t{species_base}\t"
                f"{','.join(p.resistance_bases)}\t{p.drug}\t{p.evidence_tier}\n"
            )
            if species_base != p.ref_base.upper():
                print(
                    f"  note {organism} E.coli {p.ecoli_position}: "
                    f"species base {species_base} != E. coli {p.ref_base} "
                    f"(published as {p.published_as or '-'})"
                )
        bed_path.write_text(HEADER + "".join(rows))
        msg = f"{organism}: {len(rows)} targets -> {bed_path.name}"
        if skipped:
            msg += f"  (skipped, no alignment: {skipped})"
        print(msg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
