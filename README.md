# linezolid-amr

**Integrated AMR profiling and 23S rRNA linezolid-resistance frequency analysis.**

`linezolid-amr` runs two analyses from one command and combines them into a single report:

1. **AMRFinderPlus** on an assembly — detects all AMR genes and known point mutations (cfr, optrA, poxtA, ribosomal protein L3/L4/L22 mutations, ...).
2. **23S rRNA mutation-frequency** on raw paired reads — maps to a species-specific 23S reference and reports allele frequencies at the canonical linezolid-resistance positions (G2576T, G2505A, G2447T, C2534T, A2503G, ...). This catches **heteroresistance**: mutations present in only some of the multiple 23S rRNA operon copies, which assembly-based callers miss.

The same `--organism` flag drives both AMRFinderPlus species mode and 23S reference selection. Mutations are reported in **E. coli K-12 23S numbering** (the clinical-literature convention) via a pairwise-alignment-derived position map.

## Supported organisms (23S analysis)

| Organism                       | 23S copies | Reference (1st copy)                       |
|--------------------------------|:----------:|--------------------------------------------|
| Staphylococcus aureus          | 5          | NCTC 8325 (NC_007795.1) — 2,923 bp         |
| Enterococcus faecalis          | 4          | V583 (NC_004668.1) — 2,913 bp              |
| Enterococcus faecium           | 6          | DO (NC_017960.1) — 2,914 bp                |
| Streptococcus pneumoniae       | 4          | TIGR4 (NC_003028.3) — 2,902 bp             |

AMRFinderPlus supports many more organisms; the 23S step will simply be skipped for organisms outside this list.

## Install

### From bioconda (once the recipe is accepted)

```bash
conda install -n linezolid-amr -c bioconda -c conda-forge linezolid-amr
conda activate linezolid-amr
```

### From source

```bash
git clone https://github.com/iowa69/linezolid-amr
cd linezolid-amr
conda env create -n linezolid-amr -c bioconda -c conda-forge \
    python=3.11 ncbi-amrfinderplus minimap2 samtools bcftools
conda activate linezolid-amr
pip install -e .
```

### One-time setup

```bash
# Download AMRFinderPlus database (~150 MB)
amrfinder -u

# Fetch 23S references for the species you'll analyze (~25 KB total per species)
linezolid-amr fetch-references --organism Staphylococcus_aureus
# or pull all five at once:
linezolid-amr fetch-references
```

References are cached under `$XDG_DATA_HOME/linezolid-amr/references` (defaults to `~/.local/share/linezolid-amr/references`). Override with `$LINEZOLID_AMR_REFDIR`.

## Usage

```bash
linezolid-amr run \
    --assembly  sample.fasta \
    --r1        sample_R1.fastq.gz \
    --r2        sample_R2.fastq.gz \
    --organism  Staphylococcus_aureus \
    --outdir    results/sample \
    --sample    sample \
    --threads   8
```

Subcommands:

| Command            | Purpose                                                          |
|--------------------|------------------------------------------------------------------|
| `run`              | Full pipeline (AMRFinderPlus + 23S analysis + combined report)   |
| `amrfinder`        | AMRFinderPlus step only                                          |
| `rrna23s`          | 23S analysis only (map + pileup + VCF)                           |
| `fetch-references` | Download/build references from NCBI                              |
| `list-organisms`   | List 23S-supported organisms                                     |

Important options:

* `--min-af 0.01` — minimum alt-allele frequency to report (default 1 %).
* `--min-depth 20` — minimum read depth to call a resistance allele.
* `--skip-amrfinder` / `--skip-rrna23s` — skip one of the two steps.

## Output structure

```
results/sample/
├── amrfinder/
│   ├── amrfinder.tsv                  # AMRFinderPlus raw TSV (all hits)
│   └── amrfinder.log
├── rrna23s/
│   ├── sample.23S.bam                 # sorted, indexed alignment
│   ├── sample.23S.bam.bai
│   ├── sample.23S.vcf.gz              # all variants in the 23S locus (bcftools, ploidy=1)
│   ├── sample.23S.vcf.gz.csi
│   ├── sample.23S_lzd_pileup.tsv      # allele counts/frequencies at canonical LZD positions
│   ├── sample.23S.map.log
│   └── sample.23S.vcf.log
├── sample.linezolid_amr.json          # full machine-readable report
└── sample.linezolid_amr.txt           # human-readable summary
```

The pileup TSV columns:

| column             | meaning                                                              |
|--------------------|----------------------------------------------------------------------|
| `species_position` | 1-based position in the species 23S reference                        |
| `ecoli_position`   | corresponding E. coli K-12 23S position (literature standard)        |
| `ref_base`         | wild-type base                                                       |
| `depth`            | total reads at this position                                         |
| `counts`           | per-base counts (e.g. `A=2;C=1;G=820`)                               |
| `alt_alleles`      | non-reference alleles passing thresholds (`base:count:af[*=resistance]`) |
| `is_resistance`    | `True` if a known resistance allele was observed above thresholds    |
| `drug`             | resistance phenotype the position is associated with                 |

## Interpreting heteroresistance

Because organisms like *S. aureus* and *Enterococcus* spp. have 4-6 copies of the 23S rRNA operon, low (1-50 %) alt-allele frequencies at canonical positions are biologically meaningful: they reflect the fraction of rRNA operons that carry the resistance mutation. A G2576T at 25 % AF means ~1 of 4-5 copies is mutated — a classic heteroresistant phenotype that is often missed by assembly-only callers (where the consensus base remains G).

The `summary.linezolid_resistance_call` flag in the JSON is `True` if either:

* AMRFinderPlus reports a linezolid-relevant hit (cfr / optrA / poxtA / 23S point mutation / L3-L4-L22 mutation), **or**
* The 23S pileup observes a canonical resistance allele at AF ≥ `--min-af` and depth ≥ `--min-depth`.

## Building the bioconda recipe locally

```bash
conda install -c conda-forge boa conda-build
conda-build conda-recipe/
```

## License

MIT — see [LICENSE](LICENSE).
