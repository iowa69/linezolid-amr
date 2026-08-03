"""Integrity checks on the bundled 23S references and coordinate maps.

A coordinate error here is the worst failure this tool can have: it would pile
up reads at the wrong base and report allele fractions for a position nobody
asked about, silently and with full confidence. These tests re-verify the
invariants from first principles against the shipped files.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from linezolid_amr import references as ref_mod

ORGANISMS = ref_mod.list_organisms()

# Canonical bases in E. coli K-12 rrlB. Independently anchored by classic
# landmarks (A2058/A2059 macrolide, A2451/C2452 catalytic, A2602 PTC).
ECOLI_LANDMARKS = {
    1067: "A", 2058: "A", 2059: "A", 2451: "A", 2452: "C",
    2506: "T", 2585: "T", 2602: "A", 2611: "C", 2661: "G",
}


def _read_fasta(path: Path) -> str:
    lines = path.read_text().splitlines()
    return "".join(l.strip() for l in lines[1:] if not l.startswith(">")).upper()


def _read_bed(path: Path) -> list[dict]:
    out = []
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split("\t")
        out.append({
            "chrom": f[0], "start": int(f[1]), "end": int(f[2]),
            "name": f[3], "ref_base": f[4], "resistance": f[5], "drug": f[6],
        })
    return out


@pytest.fixture(scope="module")
def ecoli_seq() -> str:
    return _read_fasta(ref_mod.ecoli_fasta_path())


# --------------------------------------------------------------------------
# The E. coli coordinate master
# --------------------------------------------------------------------------

def test_ecoli_reference_length(ecoli_seq):
    expected = ref_mod.get_ecoli_reference().expected_length_bp
    assert len(ecoli_seq) == expected


@pytest.mark.parametrize("position,base", sorted(ECOLI_LANDMARKS.items()))
def test_ecoli_landmarks(ecoli_seq, position, base):
    """Well-established 23S positions must sit where the literature puts them.

    If these drift, the reference has been re-sliced with a different offset and
    every linezolid coordinate downstream is wrong.
    """
    assert ecoli_seq[position - 1] == base


def test_loci_ref_bases_match_the_ecoli_reference(ecoli_seq):
    """Every declared E. coli ref_base must match the shipped E. coli sequence."""
    mismatches = []
    for p in ref_mod.get_linezolid_positions():
        actual = ecoli_seq[p.ecoli_position - 1]
        if actual != p.ref_base.upper():
            mismatches.append(
                f"E. coli {p.ecoli_position}: loci.json says {p.ref_base}, reference has {actual}"
            )
    assert not mismatches, "loci.json disagrees with the E. coli reference:\n  " + "\n  ".join(mismatches)


# --------------------------------------------------------------------------
# Per-species references and BED targets
# --------------------------------------------------------------------------

@pytest.mark.parametrize("organism", ORGANISMS)
def test_bed_ref_base_matches_species_fasta(organism):
    """The BED ref_base must be the base actually present at that species coordinate.

    The pileup treats every other base at this position as an alternate allele,
    so a wrong ref_base turns the species' own sequence into a spurious variant.
    """
    seq = _read_fasta(ref_mod.organism_fasta_path(organism))
    mismatches = []
    for t in _read_bed(ref_mod.organism_bed_path(organism)):
        actual = seq[t["end"] - 1]  # BED end is the 1-based position
        if actual != t["ref_base"].upper():
            mismatches.append(
                f"{t['name']} at species {t['end']}: BED says {t['ref_base']}, FASTA has {actual}"
            )
    assert not mismatches, f"{organism} BED/FASTA mismatch:\n  " + "\n  ".join(mismatches)


@pytest.mark.parametrize("organism", ORGANISMS)
def test_bed_intervals_are_single_bases(organism):
    for t in _read_bed(ref_mod.organism_bed_path(organism)):
        assert t["end"] - t["start"] == 1, f"{t['name']} is not a 1 bp interval"


@pytest.mark.parametrize("organism", ORGANISMS)
def test_bed_contig_matches_fasta_header(organism):
    header = ref_mod.organism_fasta_path(organism).read_text().splitlines()[0][1:].split()[0]
    for t in _read_bed(ref_mod.organism_bed_path(organism)):
        assert t["chrom"] == header, f"BED contig {t['chrom']} != FASTA header {header}"


@pytest.mark.parametrize("organism", ORGANISMS)
def test_bed_covers_every_declared_position(organism):
    declared = {p.ecoli_position for p in ref_mod.get_linezolid_positions()}
    bed_positions = set()
    for t in _read_bed(ref_mod.organism_bed_path(organism)):
        # name looks like 23S_E2576_G_to_T/U/A/C
        bed_positions.add(int(t["name"].split("_")[1][1:]))
    assert bed_positions == declared


# --------------------------------------------------------------------------
# Position maps
# --------------------------------------------------------------------------

@pytest.mark.parametrize("organism", ORGANISMS)
def test_position_map_is_internally_consistent(organism, ecoli_seq):
    """Both base columns in the map must agree with the two FASTAs it links."""
    species_seq = _read_fasta(ref_mod.organism_fasta_path(organism))
    path = ref_mod.organism_position_map_path(organism)
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    assert header == ["ecoli_position", "species_position", "species_base", "ecoli_base"]

    ecoli_errors = species_errors = 0
    seen_ecoli = []
    for line in lines[1:]:
        e_pos, sp_pos, sp_base, e_base = line.split("\t")
        e_pos = int(e_pos)
        seen_ecoli.append(e_pos)
        if ecoli_seq[e_pos - 1] != e_base:
            ecoli_errors += 1
        if sp_pos not in (".", ""):
            if species_seq[int(sp_pos) - 1] != sp_base:
                species_errors += 1
    assert ecoli_errors == 0, f"{organism}: {ecoli_errors} rows disagree with the E. coli FASTA"
    assert species_errors == 0, f"{organism}: {species_errors} rows disagree with the species FASTA"
    assert seen_ecoli == list(range(1, len(ecoli_seq) + 1)), "map must cover E. coli 1..N exactly once"


@pytest.mark.parametrize("organism", ORGANISMS)
def test_bed_positions_round_trip_through_the_map(organism):
    """The species coordinate in the BED must be what the map assigns to that E. coli position."""
    path = ref_mod.organism_position_map_path(organism)
    ecoli_to_species: dict[int, str] = {}
    for line in path.read_text().splitlines()[1:]:
        e_pos, sp_pos, _sp_base, _e_base = line.split("\t")
        ecoli_to_species[int(e_pos)] = sp_pos

    for t in _read_bed(ref_mod.organism_bed_path(organism)):
        ecoli_pos = int(t["name"].split("_")[1][1:])
        expected = ecoli_to_species[ecoli_pos]
        assert expected == str(t["end"]), (
            f"{organism} {t['name']}: BED says species {t['end']}, map says {expected}"
        )


@pytest.mark.parametrize("organism", ORGANISMS)
def test_rrn_copy_number_is_declared(organism):
    """Copy number sets the operon grid used to interpret allele fractions."""
    n = ref_mod.rrn_copy_number(organism)
    assert 2 <= n <= 10, f"{organism}: implausible rrn copy number {n}"


# --------------------------------------------------------------------------
# Curation: evidence tiers and literature notation
# --------------------------------------------------------------------------

def test_every_position_has_a_known_tier():
    for p in ref_mod.get_linezolid_positions():
        assert p.evidence_tier in ref_mod.EVIDENCE_TIERS, (
            f"E. coli {p.ecoli_position} has unknown tier {p.evidence_tier!r}"
        )


def test_every_position_has_a_pmid():
    for p in ref_mod.get_linezolid_positions():
        assert p.pmid, f"position {p.ecoli_position} has no PMID citation"
        for pmid in p.pmid:
            assert pmid.isdigit() and len(pmid) >= 7, f"invalid PMID {pmid}"


def test_g2576_is_the_established_anchor():
    """The dominant clinical mutation must be present and callable by default."""
    g2576 = [p for p in ref_mod.get_linezolid_positions() if p.ecoli_position == 2576]
    assert len(g2576) == 1
    p = g2576[0]
    assert p.ref_base == "G"
    assert "T" in p.resistance_bases
    assert p.evidence_tier == ref_mod.TIER_ESTABLISHED
    assert ref_mod.tier_is_callable(p.evidence_tier)


def test_a2572_screens_for_the_reported_substitution():
    """Long & Vester Table 1 lists A2572U. Screening for A2572G was a curation bug."""
    p = next(x for x in ref_mod.get_linezolid_positions() if x.ecoli_position == 2572)
    assert "T" in p.resistance_bases or "U" in p.resistance_bases
    assert "G" not in p.resistance_bases


def test_hybrid_notation_positions_carry_a_caveat():
    """2534 and 2571 are published using a non-E. coli wild-type letter.

    Both must record the discrepancy explicitly, or a reader comparing the
    output against the literature will conclude the tool is miscalling.
    """
    for pos in (2534, 2571):
        p = next(x for x in ref_mod.get_linezolid_positions() if x.ecoli_position == pos)
        assert p.numbering_caveat, f"position {pos} needs a numbering caveat"
        # The published letter must differ from the true E. coli base.
        assert p.published_as
        assert p.published_as[0] != p.ref_base


def test_experimental_positions_do_not_call_by_default():
    """Lab-only substitutions are reported but must not manufacture a positive call."""
    experimental = [
        p for p in ref_mod.get_linezolid_positions()
        if p.evidence_tier == ref_mod.TIER_EXPERIMENTAL
    ]
    assert experimental, "expected some experimental-tier positions"
    for p in experimental:
        assert not ref_mod.tier_is_callable(p.evidence_tier)


def test_tier_filter_narrows_the_position_set():
    established = ref_mod.get_linezolid_positions(min_tier=ref_mod.TIER_ESTABLISHED)
    default = ref_mod.get_linezolid_positions(min_tier=ref_mod.DEFAULT_EVIDENCE_TIER)
    everything = ref_mod.get_linezolid_positions()
    assert len(established) < len(default) < len(everything)


def test_tier_ordering():
    assert ref_mod.tier_is_callable("established", "established")
    assert not ref_mod.tier_is_callable("associated", "established")
    assert ref_mod.tier_is_callable("associated", "associated")
    assert ref_mod.tier_is_callable("experimental", "experimental")
    assert not ref_mod.tier_is_callable("experimental", "associated")


@pytest.mark.parametrize("organism", ORGANISMS)
def test_bed_carries_the_evidence_tier(organism):
    for line in ref_mod.organism_bed_path(organism).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        fields = line.split("\t")
        assert len(fields) >= 8, f"BED row missing the tier column: {line}"
        assert fields[7] in ref_mod.EVIDENCE_TIERS


def test_excluded_positions_are_documented():
    """Deliberate exclusions are recorded so a reviewer can see they were considered."""
    data = ref_mod.load_loci()
    assert data.get("excluded_positions"), "loci.json should record what was excluded and why"
    for item in data["excluded_positions"]:
        assert item.get("reason")


def test_staging_copies_reference_out_of_the_package(tmp_path):
    organism = ORGANISMS[0]
    staged = ref_mod.stage_reference(organism, tmp_path)
    src = ref_mod.organism_fasta_path(organism)
    assert staged.parent == tmp_path
    assert staged.read_bytes() == src.read_bytes()


def test_staging_is_idempotent(tmp_path):
    organism = ORGANISMS[0]
    first = ref_mod.stage_reference(organism, tmp_path)
    first.with_suffix(first.suffix + ".fai").write_text("stale\n")
    second = ref_mod.stage_reference(organism, tmp_path)
    assert first == second
    assert second.exists()
