"""linezolid-amr — command-line interface."""

from __future__ import annotations

import os
import sys
from pathlib import Path

import click

from linezolid_amr import __version__
from linezolid_amr import amrfinder as amr_mod
from linezolid_amr import fetch_references as fetch_mod
from linezolid_amr import internal_mlst as mlst_mod
from linezolid_amr import references as ref_mod
from linezolid_amr import reporting as report_mod
from linezolid_amr import rrna23s as rrna_mod
from linezolid_amr import summary as summary_mod


DEFAULT_THREADS = max(os.cpu_count() or 4, 1)

_SUPPORTED_ORGS = ref_mod.list_organisms()

ORGANISMS_HELP = "Supported organisms (23S + MLST): " + ", ".join(_SUPPORTED_ORGS) + "."

MAIN_DESCRIPTION = f"""linezolid-amr — integrated AMR profiling + 23S rRNA linezolid-resistance analysis.

\b
SUPPORTED ORGANISMS (both 23S heteroresistance and MLST sequence typing):
  • Staphylococcus aureus
  • Enterococcus faecalis
  • Enterococcus faecium
  • Streptococcus pneumoniae

\b
PIPELINE STAGES:
  1. MLST → infers organism + Sequence Type from the assembly (in-house, no external mlst tool)
  2. AMRFinderPlus → all AMR / virulence / stress genes + canonical point mutations
  3. 23S rRNA read pileup → allele frequencies at canonical linezolid-resistance positions
                          (catches heteroresistance via multi-operon mixed reads)

For organisms outside the four above the 23S and MLST steps are skipped and
only AMRFinderPlus runs (organism passed via -O / --organism)."""


_BANNER = r"""
   _ _                           _ _     _
  | (_)                         | (_)   | |
  | |_ _ __   ___ _______  _ __ | |_  __| |   __ _ _ __ ___  _ __
  | | | '_ \ / _ \_  / _ \| '_ \| | |/ _` |  / _` | '_ ` _ \| '__|
  | | | | | |  __// / (_) | | | | | | (_| | | (_| | | | | | | |
  |_|_|_| |_|\___/___\___/|_| |_|_|_|\__,_|  \__,_|_| |_| |_|_|
"""


def _print_banner() -> None:
    """Show banner + supported organism list before any run."""
    click.echo(_BANNER)
    click.echo("Supported organisms: " + ", ".join(_SUPPORTED_ORGS))
    click.echo("")


# ---------------- helpers for folder mode ---------------- #

_ASSEMBLY_EXTS = (".fasta", ".fa", ".fas", ".fna")
_FASTQ_EXTS    = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
_R1_PATTERNS   = ("_R1_001", "_R1", "_1")
_R2_PATTERNS   = ("_R2_001", "_R2", "_2")


def _strip_assembly_ext(name: str) -> str | None:
    for ext in _ASSEMBLY_EXTS:
        if name.lower().endswith(ext):
            return name[: -len(ext)]
    return None


def _find_paired_reads(basename: str, folder: Path) -> tuple[Path, Path | None] | None:
    candidates = []
    for f in folder.iterdir():
        if not f.is_file():
            continue
        n = f.name
        if not any(n.lower().endswith(e) for e in _FASTQ_EXTS):
            continue
        candidates.append(f)

    def match(pat_set, candidates_):
        for f in candidates_:
            n = f.name
            for ext in _FASTQ_EXTS:
                if not n.lower().endswith(ext):
                    continue
                stem = n[: -len(ext)]
                for pat in pat_set:
                    if stem == basename + pat:
                        return f
        return None

    r1 = match(_R1_PATTERNS, candidates)
    if not r1:
        return None
    r2 = match(_R2_PATTERNS, candidates)
    return r1, r2


def _strip_fastq_ext(name: str) -> str | None:
    for ext in _FASTQ_EXTS:
        if name.lower().endswith(ext):
            return name[: -len(ext)]
    return None


def _read_basename(stem: str) -> str | None:
    """Return the sample basename from an R1 read filename stem.

    Strips one of the recognised R1 suffixes (``_R1_001``/``_R1``/``_1``).
    Returns None if stem is not an R1 read.
    """
    for suf in _R1_PATTERNS:
        if stem.endswith(suf):
            return stem[: -len(suf)]
    return None


def discover_samples(folder: Path, reads_only: bool = False) -> list[dict]:
    """Discover sample triples in *folder*.

    Default behaviour (reads_only=False):
      - First match assemblies (``*.fasta/.fa/.fas/.fna``), then look for paired
        FASTQs sharing the basename. Samples without paired reads are tagged
        ``missing_reads=True``.
      - Additionally, any R1/R2 pair that has *no* matching assembly is also
        emitted, but with ``assembly=None`` so the caller can decide whether to
        run reads-only mode (requires --organism).

    reads_only=True: only enumerate R1/R2 pairs; no assemblies are considered.
    """
    folder = folder.resolve()
    out: list[dict] = []
    consumed_reads: set[Path] = set()

    if not reads_only:
        for f in sorted(folder.iterdir()):
            if not f.is_file():
                continue
            base = _strip_assembly_ext(f.name)
            if base is None:
                continue
            reads = _find_paired_reads(base, folder)
            if reads is None:
                out.append({"sample": base, "assembly": f, "r1": None, "r2": None, "missing_reads": True})
                continue
            r1, r2 = reads
            consumed_reads.add(r1)
            if r2:
                consumed_reads.add(r2)
            out.append({"sample": base, "assembly": f, "r1": r1, "r2": r2})

    # Reads-only samples: every R1 file whose pair hasn't been claimed by an assembly
    seen_bases = {s["sample"] for s in out}
    for f in sorted(folder.iterdir()):
        if not f.is_file() or f in consumed_reads:
            continue
        stem = _strip_fastq_ext(f.name)
        if stem is None:
            continue
        base = _read_basename(stem)
        if base is None or base in seen_bases:
            continue
        reads = _find_paired_reads(base, folder)
        if reads is None:
            continue
        r1, r2 = reads
        if r1 in consumed_reads:
            continue
        out.append({"sample": base, "assembly": None, "r1": r1, "r2": r2, "reads_only": True})
        seen_bases.add(base)
    return out


# ---------------- shared pipeline ---------------- #

def _fmt_alleles(res) -> str:
    """Seemann-style allele string ``arcC(3)|aroE(35)|...`` ordered by scheme loci."""
    if res is None:
        return ""
    return "|".join(f"{l}({res.alleles[l]})" for l in res.alleles)


def _resolve_organism(
    user_organism: str | None,
    assembly: Path,
    threads: int,
) -> tuple[str | None, str | None, str | None, str]:
    """Return (organism, mlst_scheme, st, mlst_alleles_string). In-house MLST infers organism."""
    mlst_scheme = None
    st = None
    mlst_organism = None
    mlst_alleles = ""

    if mlst_mod.blastn_available() and mlst_mod.all_schemes_available():
        try:
            mlst_organism, res = mlst_mod.infer_organism(assembly, threads=threads)
            mlst_scheme, st = res.scheme, res.st
            mlst_alleles = _fmt_alleles(res)
        except mlst_mod.MlstNoMatch as e:
            if not user_organism:
                raise click.ClickException(str(e))
            click.echo(f"   MLST: no scheme matched; using user --organism", err=True)
        except mlst_mod.MlstUnsupportedOrganism as e:
            if not user_organism:
                raise click.ClickException(str(e))
            click.echo(f"   MLST: {e} (using user --organism)", err=True)
    elif not user_organism:
        msg = (
            "Internal MLST cannot run: "
            + ("`blastn` not on PATH (install via `conda install -c bioconda blast`). " if not mlst_mod.blastn_available() else "")
            + ("MLST schemes not bundled in this build — run `linezolid-amr fetch-mlst-schemes`. " if not mlst_mod.all_schemes_available() else "")
            + "Pass --organism manually or fix the above."
        )
        raise click.ClickException(msg)

    organism = user_organism or mlst_organism
    if user_organism and mlst_organism and user_organism != mlst_organism:
        click.echo(
            f"   ⚠ MLST suggests {mlst_organism} but --organism is {user_organism}. Using --organism.",
            err=True,
        )
    return organism, mlst_scheme, st, mlst_alleles


def _run_single(
    sample: str, assembly: Path | None, r1: Path, r2: Path | None,
    user_organism: str | None, outdir: Path, threads: int,
    min_af: float, min_depth: int, plus: bool,
    skip_amrfinder: bool, skip_rrna23s: bool,
) -> dict:
    """Per-sample pipeline.

    Reads-only mode: when ``assembly is None`` the MLST and AMRFinderPlus
    stages are skipped; only the 23S read-level analysis runs. The user must
    supply ``--organism`` in that case.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    click.echo(f"\n=== {sample} ===")

    mlst_scheme = mlst_alleles = None
    st = None
    if assembly is not None:
        click.echo(">> MLST / organism inference...")
        organism, mlst_scheme, st, mlst_alleles = _resolve_organism(user_organism, assembly, threads)
        click.echo(f"   organism: {organism}   MLST scheme: {mlst_scheme or '-'}   ST: {st or '-'}")
        if mlst_alleles:
            click.echo(f"   alleles: {mlst_alleles}")
    else:
        # 23S-only mode: no assembly → no MLST, AMRFinderPlus skipped, --organism required.
        if not user_organism:
            raise click.ClickException(
                f"Reads-only mode requires --organism (no assembly supplied). "
                f"Pass one of: {', '.join(_SUPPORTED_ORGS)}."
            )
        organism = user_organism
        click.echo(f">> Reads-only mode (no assembly): organism={organism}")
        skip_amrfinder = True  # nothing to scan

    parameters = {
        "sample": sample, "assembly": str(assembly) if assembly else None,
        "r1": str(r1), "r2": str(r2) if r2 else None,
        "organism": organism, "mlst_scheme": mlst_scheme, "st": st,
        "threads": threads, "min_af": min_af, "min_depth": min_depth, "plus": plus,
        "reads_only": assembly is None,
    }

    amr_hits = []
    if not skip_amrfinder and assembly is not None:
        click.echo(">> Running AMRFinderPlus...")
        tsv = amr_mod.run_amrfinder(assembly, organism, outdir / "amrfinder", threads=threads, plus=plus)
        amr_hits = amr_mod.parse_amrfinder_tsv(tsv)
        click.echo(f"   {len(amr_hits)} hits, {len(amr_mod.linezolid_relevant_hits(amr_hits))} linezolid-relevant")

    pileup_calls = []
    vcf = None
    if not skip_rrna23s:
        if organism not in _SUPPORTED_ORGS:
            click.echo(f"   organism '{organism}' not supported for 23S analysis — skipping", err=True)
        else:
            click.echo(">> Running 23S rRNA analysis...")
            rrna_outdir = outdir / "rrna23s"
            fasta = ref_mod.organism_fasta_path(organism)
            ref_mod.ensure_references_available(organism)
            bam = rrna_mod.map_reads(fasta, r1, r2, rrna_outdir, threads=threads, sample=sample)
            pileup_calls = rrna_mod.pileup_at_positions(bam, fasta, organism, min_af=min_af, min_depth=min_depth)
            rrna_mod.write_pileup_tsv(pileup_calls, rrna_outdir / f"{sample}.23S_lzd_pileup.tsv")
            vcf = rrna_mod.call_full_vcf(bam, fasta, rrna_outdir, sample=sample, threads=threads)
            click.echo(f"   {len(pileup_calls)} positions; {sum(1 for c in pileup_calls if c.is_resistance)} with resistance allele")

    report = report_mod.build_report(
        sample=sample, organism=organism, amr_hits=amr_hits,
        pileup_calls=pileup_calls, vcf_path=vcf, parameters=parameters,
    )
    report["mlst"] = {"scheme": mlst_scheme, "ST": st}

    json_path = outdir / f"{sample}.linezolid_amr.json"
    txt_path = outdir / f"{sample}.linezolid_amr.txt"
    report_mod.write_json(report, json_path)
    report_mod.write_text_summary(report, txt_path)

    lzd_call = bool(report["summary"]["linezolid_resistance_call"])
    wide_row = summary_mod.build_wide_row(
        sample, organism, st, mlst_scheme, mlst_alleles,
        amr_hits, pileup_calls, lzd_call,
    )
    long_rows = summary_mod.build_long_rows(sample, organism, st, mlst_scheme, amr_hits, pileup_calls)
    summary_mod.write_wide_csv([wide_row], outdir / f"{sample}.summary_wide.csv")
    summary_mod.write_long_csv(long_rows, outdir / f"{sample}.summary_long.csv")

    return {
        "sample": sample, "organism": organism, "mlst_scheme": mlst_scheme, "ST": st,
        "linezolid_call": lzd_call, "wide_row": wide_row, "long_rows": long_rows,
        "json": json_path, "txt": txt_path,
    }


# ---------------- click commands ---------------- #

@click.group(
    context_settings={"help_option_names": ["-h", "--help"], "max_content_width": 100},
    help=MAIN_DESCRIPTION,
)
@click.version_option(__version__, prog_name="linezolid-amr")
def main() -> None:
    pass


@main.command("list-organisms", help="List organisms supported by the 23S + MLST steps.")
def list_organisms_cmd() -> None:
    for org in _SUPPORTED_ORGS:
        try:
            o = ref_mod.get_organism(org)
            click.echo(f"{org}\t{o.description}")
        except KeyError:
            click.echo(org)


@main.command("fetch-references",
              help=f"Download per-species 23S references and build E. coli position maps. "
                   f"{ORGANISMS_HELP}")
@click.option("-O", "--organism", "organisms", multiple=True,
              help="Limit to one or more organisms (default: all supported).")
@click.option("--api-key", envvar="NCBI_API_KEY", default=None,
              help="Optional NCBI E-utilities API key.")
@click.option("--force", is_flag=True, help="Re-download even if cached.")
def fetch_references_cmd(organisms, api_key, force):
    click.echo(f"Cache directory: {ref_mod.cache_dir()}")
    results = fetch_mod.fetch_all(
        organisms=list(organisms) if organisms else None,
        api_key=api_key, force=force,
    )
    failures = 0
    for r in results:
        if r.get("status") == "error":
            failures += 1
            click.echo(f"  [ERR ] {r['organism']}: {r.get('error')}", err=True)
        else:
            click.echo(f"  [{r['status'].upper():>6}] {r['organism']}")
    if failures:
        sys.exit(1)


@main.command("fetch-mlst-schemes",
              help=f"Download/refresh PubMLST schemes for the supported organisms. "
                   f"Schemes are bundled in the package by default — use this only if you "
                   f"want fresher allele sets than the release tag carries. {ORGANISMS_HELP}")
@click.option("-O", "--organism", "organisms", multiple=True,
              help="Limit to specific organisms (default: all supported).")
@click.option("--outdir", type=click.Path(path_type=Path), default=None,
              help="Output dir (default: $LINEZOLID_AMR_MLST_DIR or ~/.local/share/linezolid-amr/mlst_schemes/).")
@click.option("--force", is_flag=True, help="Re-download even if cached.")
def fetch_mlst_cmd(organisms, outdir, force):
    targets = list(organisms) if organisms else list(mlst_mod.PUBMLST_DBS)
    click.echo(f"Fetching PubMLST schemes ({len(targets)} organism(s))...")
    results = mlst_mod.fetch_pubmlst_schemes(organisms=targets, out_root=outdir, force=force)
    failures = 0
    for r in results:
        if r.get("status") == "error":
            failures += 1
            click.echo(f"  [ERR ] {r['organism']}: {r.get('error')}", err=True)
        elif r.get("status") == "fetched":
            click.echo(f"  [FETCHED] {r['organism']} → {r['n_alleles']} alleles in {r['dir']}")
        else:
            click.echo(f"  [{r['status'].upper():>6}] {r['organism']}")
    if failures:
        sys.exit(1)


# Shared option set
def _shared_run_options(f):
    f = click.option("-O", "--organism", default=None,
                     help=f"AMRFinderPlus organism. If omitted, inferred via in-house MLST. {ORGANISMS_HELP}")(f)
    f = click.option("-t", "--threads", default=DEFAULT_THREADS, show_default=True, type=int,
                     help="Threads (default: all available CPUs).")(f)
    f = click.option("--min-af", default=rrna_mod.DEFAULT_MIN_AF, show_default=True, type=float,
                     help="Minimum alt-allele frequency to call a 23S resistance allele.")(f)
    f = click.option("--min-depth", default=rrna_mod.DEFAULT_MIN_DEPTH, show_default=True, type=int,
                     help="Minimum read depth at a 23S position.")(f)
    f = click.option("--plus", is_flag=True, default=False,
                     help="Pass AMRFinderPlus --plus (stress/virulence/biocide search).")(f)
    f = click.option("--skip-amrfinder", is_flag=True, help="Skip AMRFinderPlus step.")(f)
    f = click.option("--skip-rrna23s", is_flag=True, help="Skip 23S analysis step.")(f)
    return f


@main.command("run",
              help=f"Run the pipeline on a single sample.\n\n"
                   f"With -a/--assembly: full pipeline (MLST + AMRFinderPlus + 23S).\n"
                   f"Without -a: reads-only 23S linezolid-resistance call; -O/--organism required.\n\n"
                   f"{ORGANISMS_HELP}")
@click.option("-a", "--assembly", required=False, default=None,
              type=click.Path(exists=True, path_type=Path),
              help="Genome assembly FASTA. Omit for reads-only 23S analysis (--organism required).")
@click.option("-1", "--r1", required=True, type=click.Path(exists=True, path_type=Path),
              help="Forward reads (FASTQ, optionally gzip).")
@click.option("-2", "--r2", default=None, type=click.Path(exists=True, path_type=Path),
              help="Reverse reads (optional).")
@click.option("-o", "--outdir", required=True, type=click.Path(path_type=Path),
              help="Output directory.")
@click.option("-s", "--sample", default=None,
              help="Sample name (default: assembly basename or R1 basename).")
@_shared_run_options
def run_cmd(assembly, r1, r2, outdir, sample, organism, threads, min_af, min_depth, plus,
            skip_amrfinder, skip_rrna23s):
    _print_banner()
    if sample is None:
        if assembly is not None:
            sample = assembly.stem
        else:
            # derive from R1 basename, stripping common suffixes
            base = r1.name
            for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
                if base.lower().endswith(ext):
                    base = base[: -len(ext)]
                    break
            for suf in ("_R1_001", "_R1", "_1"):
                if base.endswith(suf):
                    base = base[: -len(suf)]
                    break
            sample = base
    res = _run_single(
        sample, assembly, r1, r2, organism, outdir, threads, min_af, min_depth, plus,
        skip_amrfinder, skip_rrna23s,
    )
    click.echo("")
    click.echo(f"Report:        {res['json']}")
    click.echo(f"Summary (txt): {res['txt']}")
    click.echo(f"Summary CSVs:  {outdir}/{sample}.summary_{{wide,long}}.csv")
    click.echo(f"Linezolid resistance call: {'POSITIVE' if res['linezolid_call'] else 'negative'}")


@main.command("folder",
              help=f"Run the pipeline on every sample in a folder, emit a cohort summary.\n\n"
                   f"Discovery: *.fasta/*.fa/*.fas/*.fna paired with FASTQs by basename "
                   f"(R1/R2/_1/_2/_001 × .fastq[.gz]/.fq[.gz]).\n\n"
                   f"With --reads-only: ignore assemblies entirely, batch every R1/R2 pair "
                   f"in reads-only 23S mode (--organism required).\n\n{ORGANISMS_HELP}")
@click.option("-i", "--input", "input_dir", required=True,
              type=click.Path(exists=True, file_okay=False, path_type=Path),
              help="Folder containing assemblies + paired FASTQs.")
@click.option("-o", "--outdir", required=True, type=click.Path(path_type=Path),
              help="Output directory; one subfolder per sample.")
@click.option("--reads-only", is_flag=True, default=False,
              help="Skip assembly discovery; process every R1/R2 pair in 23S-only mode (requires --organism).")
@_shared_run_options
def folder_cmd(input_dir, outdir, organism, threads, min_af, min_depth, plus,
               skip_amrfinder, skip_rrna23s, reads_only):
    _print_banner()
    samples = discover_samples(input_dir, reads_only=reads_only)
    if not samples:
        raise click.ClickException(
            f"No samples found in {input_dir}. "
            f"{'Looked for FASTQ R1/R2 pairs.' if reads_only else 'Looked for *.fasta + paired FASTQs and bare R1/R2 pairs.'}"
        )

    valid = [s for s in samples if not s.get("missing_reads")]
    skipped = [s for s in samples if s.get("missing_reads")]
    n_full = sum(1 for s in valid if s.get("assembly") is not None)
    n_reads_only = sum(1 for s in valid if s.get("assembly") is None)
    if n_reads_only and not organism:
        raise click.ClickException(
            f"{n_reads_only} sample(s) have only FASTQ reads (no matching assembly). "
            f"Pass -O/--organism for reads-only 23S mode or supply the assemblies."
        )
    click.echo(
        f"Discovered {len(valid)} samples "
        f"({n_full} full pipeline, {n_reads_only} reads-only; skipping {len(skipped)} without paired reads)."
    )
    for s in skipped:
        click.echo(f"   ⚠ {s['sample']}: paired reads not found in {input_dir}", err=True)

    outdir.mkdir(parents=True, exist_ok=True)
    all_wide_rows: list[dict] = []
    all_long_rows: list[dict] = []
    for s in valid:
        sample_dir = outdir / s["sample"]
        try:
            res = _run_single(
                s["sample"], s["assembly"], s["r1"], s["r2"],
                organism, sample_dir, threads, min_af, min_depth, plus,
                skip_amrfinder, skip_rrna23s,
            )
            all_wide_rows.append(res["wide_row"])
            all_long_rows.extend(res["long_rows"])
        except Exception as e:  # noqa: BLE001
            click.echo(f"   ✗ {s['sample']} failed: {e}", err=True)

    summary_mod.write_wide_csv(all_wide_rows, outdir / "ALL_samples.summary_wide.csv")
    summary_mod.write_long_csv(all_long_rows, outdir / "ALL_samples.summary_long.csv")

    click.echo("")
    click.echo(f"Cohort summary (wide): {outdir}/ALL_samples.summary_wide.csv")
    click.echo(f"Cohort summary (long): {outdir}/ALL_samples.summary_long.csv")
    n_pos = sum(1 for r in all_wide_rows if r.get("linezolid_call") == "POS")
    click.echo(f"Linezolid-positive samples: {n_pos} / {len(all_wide_rows)}")


if __name__ == "__main__":
    main()
