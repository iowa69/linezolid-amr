"""linezolid-amr command-line interface."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import click

from linezolid_amr import __version__
from linezolid_amr import amrfinder as amr_mod
from linezolid_amr import fetch_references as fetch_mod
from linezolid_amr import references as ref_mod
from linezolid_amr import reporting as report_mod
from linezolid_amr import rrna23s as rrna_mod


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.version_option(__version__, prog_name="linezolid-amr")
def main() -> None:
    """Integrated AMR profiling and 23S rRNA linezolid-resistance frequency analysis.

    Combines NCBI AMRFinderPlus (gene/mutation detection from an assembly) with
    species-aware mapping of raw reads to 23S rRNA references to detect
    heteroresistant linezolid-resistance mutations (e.g. G2576T).
    """


@main.command("list-organisms")
def list_organisms_cmd() -> None:
    """List organisms supported by the 23S rRNA analysis."""
    for org in ref_mod.list_organisms():
        try:
            o = ref_mod.get_organism(org)
            click.echo(f"{org}\t{o.description}")
        except KeyError:
            click.echo(org)


@main.command("fetch-references")
@click.option(
    "--organism",
    "organisms",
    multiple=True,
    help="Limit to one or more organisms (default: all supported).",
)
@click.option(
    "--api-key",
    envvar="NCBI_API_KEY",
    default=None,
    help="Optional NCBI E-utilities API key (also reads $NCBI_API_KEY).",
)
@click.option("--force", is_flag=True, help="Re-download even if cached.")
def fetch_references_cmd(organisms: tuple[str, ...], api_key: str | None, force: bool) -> None:
    """Download per-species 23S references and build E. coli position maps."""
    click.echo(f"Cache directory: {ref_mod.cache_dir()}")
    results = fetch_mod.fetch_all(
        organisms=list(organisms) if organisms else None,
        api_key=api_key,
        force=force,
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


@main.command("amrfinder")
@click.option("--assembly", required=True, type=click.Path(exists=True, path_type=Path))
@click.option("--organism", required=True)
@click.option("--outdir", required=True, type=click.Path(path_type=Path))
@click.option("--threads", default=4, show_default=True, type=int)
def amrfinder_cmd(assembly: Path, organism: str, outdir: Path, threads: int) -> None:
    """Run only the AMRFinderPlus analysis."""
    tsv = amr_mod.run_amrfinder(assembly, organism, outdir, threads=threads)
    hits = amr_mod.parse_amrfinder_tsv(tsv)
    click.echo(f"AMRFinderPlus: {len(hits)} hits ({len(amr_mod.linezolid_relevant_hits(hits))} linezolid-relevant) -> {tsv}")


@main.command("rrna23s")
@click.option("--r1", required=True, type=click.Path(exists=True, path_type=Path))
@click.option("--r2", default=None, type=click.Path(exists=True, path_type=Path))
@click.option("--organism", required=True)
@click.option("--outdir", required=True, type=click.Path(path_type=Path))
@click.option("--sample", default="sample", show_default=True)
@click.option("--threads", default=4, show_default=True, type=int)
@click.option("--min-af", default=rrna_mod.DEFAULT_MIN_AF, show_default=True, type=float)
@click.option("--min-depth", default=rrna_mod.DEFAULT_MIN_DEPTH, show_default=True, type=int)
def rrna23s_cmd(
    r1: Path, r2: Path | None, organism: str, outdir: Path, sample: str,
    threads: int, min_af: float, min_depth: int,
) -> None:
    """Run only the 23S rRNA frequency analysis (mapping + pileup + VCF)."""
    fasta = ref_mod.organism_fasta_path(organism)
    ref_mod.ensure_references_available(organism)
    bam = rrna_mod.map_reads(fasta, r1, r2, outdir, threads=threads, sample=sample)
    calls = rrna_mod.pileup_at_positions(bam, fasta, organism, min_af=min_af, min_depth=min_depth)
    rrna_mod.write_pileup_tsv(calls, outdir / f"{sample}.23S_lzd_pileup.tsv")
    vcf = rrna_mod.call_full_vcf(bam, fasta, outdir, sample=sample, threads=threads)
    click.echo(f"23S pileup: {sum(1 for c in calls if c.is_resistance)} resistance positions; VCF -> {vcf}")


@main.command("run")
@click.option("--assembly", required=True, type=click.Path(exists=True, path_type=Path))
@click.option("--r1", required=True, type=click.Path(exists=True, path_type=Path))
@click.option("--r2", default=None, type=click.Path(exists=True, path_type=Path))
@click.option("--organism", required=True)
@click.option("--outdir", required=True, type=click.Path(path_type=Path))
@click.option("--sample", default="sample", show_default=True)
@click.option("--threads", default=4, show_default=True, type=int)
@click.option("--min-af", default=rrna_mod.DEFAULT_MIN_AF, show_default=True, type=float)
@click.option("--min-depth", default=rrna_mod.DEFAULT_MIN_DEPTH, show_default=True, type=int)
@click.option("--skip-amrfinder", is_flag=True, help="Skip AMRFinderPlus step.")
@click.option("--skip-rrna23s", is_flag=True, help="Skip 23S analysis step.")
def run_cmd(
    assembly: Path, r1: Path, r2: Path | None, organism: str, outdir: Path,
    sample: str, threads: int, min_af: float, min_depth: int,
    skip_amrfinder: bool, skip_rrna23s: bool,
) -> None:
    """Run the full pipeline: AMRFinderPlus on assembly + 23S analysis on reads."""
    outdir.mkdir(parents=True, exist_ok=True)
    parameters = {
        "assembly": str(assembly), "r1": str(r1), "r2": str(r2) if r2 else None,
        "organism": organism, "sample": sample, "threads": threads,
        "min_af": min_af, "min_depth": min_depth,
    }

    amr_hits = None
    if not skip_amrfinder:
        click.echo(">> Running AMRFinderPlus...")
        amr_outdir = outdir / "amrfinder"
        tsv = amr_mod.run_amrfinder(assembly, organism, amr_outdir, threads=threads)
        amr_hits = amr_mod.parse_amrfinder_tsv(tsv)
        click.echo(f"   {len(amr_hits)} hits, {len(amr_mod.linezolid_relevant_hits(amr_hits))} linezolid-relevant")

    pileup_calls = None
    vcf = None
    if not skip_rrna23s:
        if organism not in ref_mod.list_organisms():
            click.echo(
                f"   organism '{organism}' not supported for 23S analysis; skipping. "
                f"Supported: {', '.join(ref_mod.list_organisms())}",
                err=True,
            )
        else:
            click.echo(">> Running 23S rRNA analysis...")
            rrna_outdir = outdir / "rrna23s"
            fasta = ref_mod.organism_fasta_path(organism)
            ref_mod.ensure_references_available(organism)
            bam = rrna_mod.map_reads(fasta, r1, r2, rrna_outdir, threads=threads, sample=sample)
            pileup_calls = rrna_mod.pileup_at_positions(
                bam, fasta, organism, min_af=min_af, min_depth=min_depth
            )
            rrna_mod.write_pileup_tsv(pileup_calls, rrna_outdir / f"{sample}.23S_lzd_pileup.tsv")
            vcf = rrna_mod.call_full_vcf(bam, fasta, rrna_outdir, sample=sample, threads=threads)
            n_res = sum(1 for c in pileup_calls if c.is_resistance)
            click.echo(f"   {len(pileup_calls)} positions; {n_res} with resistance allele")

    click.echo(">> Building report...")
    report = report_mod.build_report(
        sample=sample, organism=organism, amr_hits=amr_hits,
        pileup_calls=pileup_calls, vcf_path=vcf, parameters=parameters,
    )
    json_path = outdir / f"{sample}.linezolid_amr.json"
    txt_path = outdir / f"{sample}.linezolid_amr.txt"
    report_mod.write_json(report, json_path)
    report_mod.write_text_summary(report, txt_path)

    click.echo("")
    click.echo(f"Report:  {json_path}")
    click.echo(f"Summary: {txt_path}")
    click.echo(
        f"Linezolid resistance call: {report['summary']['linezolid_resistance_call']}"
    )


if __name__ == "__main__":
    main()
