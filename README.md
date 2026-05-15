<div align="center">

```
   _ _                           _ _     _
  | (_)                         | (_)   | |
  | |_ _ __   ___ _______  _ __ | |_  __| |   __ _ _ __ ___  _ __
  | | | '_ \ / _ \_  / _ \| '_ \| | |/ _` |  / _` | '_ ` _ \| '__|
  | | | | | |  __// / (_) | | | | | | (_| | | (_| | | | | | | |
  |_|_|_| |_|\___/___\___/|_| |_|_|_|\__,_|  \__,_|_| |_| |_|_|
```

**🧬 Integrated AMR profiling + 23S rRNA linezolid-resistance frequency analysis**

[![CI](https://github.com/iowa69/linezolid-amr/actions/workflows/ci.yml/badge.svg)](https://github.com/iowa69/linezolid-amr/actions)
[![bioconda](https://img.shields.io/conda/dn/bioconda/linezolid-amr.svg?label=bioconda)](https://anaconda.org/bioconda/linezolid-amr)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

</div>

---

## 🎯 What it does

One command runs **three** analyses and writes a combined report:

| Step | Tool | What it gives you |
|---|---|---|
| 1️⃣ **MLST** | `mlst` (Seemann) | Sequence Type + organism auto-detection |
| 2️⃣ **AMR profile** | NCBI AMRFinderPlus | All AMR/virulence/stress genes & point mutations |
| 3️⃣ **23S heteroresistance** | minimap2 + pysam | Allele frequency at 11 canonical LZD positions |

🧪 **Supported organisms** (23S step): *Staphylococcus aureus*, *Enterococcus faecalis*, *Enterococcus faecium*, *Streptococcus pneumoniae*. Linezolid resistance is reported in **E. coli K-12 23S numbering** (clinical convention) — translated automatically.

🔬 **Why read-level?** *S. aureus* and *Enterococcus* spp. carry 4–6 copies of the 23S rRNA operon. A G2576T mutation in only a subset of copies produces heteroresistance that assembly-only callers miss. Reading allele frequencies straight off the reads is the only way to catch it.

---

## ⚡ Install

```bash
# 🐍 from bioconda (once accepted)
conda create -n linezolid-amr -c bioconda -c conda-forge linezolid-amr
conda activate linezolid-amr

# or from source
git clone https://github.com/iowa69/linezolid-amr && cd linezolid-amr
conda create -n linezolid-amr -c bioconda -c conda-forge \
    python=3.11 ncbi-amrfinderplus minimap2 samtools bcftools mlst
conda activate linezolid-amr
pip install -e .
```

One-time downloads (the tool prompts if missing):
```bash
amrfinder -u                                    # AMRFinder DB (~150 MB)
linezolid-amr fetch-references                  # 23S references (~120 KB)
```

---

## 🚀 Usage

### Single sample

```bash
linezolid-amr run -a sample.fasta -1 sample_R1.fq.gz -2 sample_R2.fq.gz \
                  -o results/sample
```

`--organism` is optional — MLST infers it. Pass `-O Enterococcus_faecium` to override.

### 🗂️ Folder mode (batch)

Drop all `*.fasta` + paired `*.fastq.gz` in one folder. Recognized read suffixes:
`_R1_001`/`_R2_001`, `_R1`/`_R2`, `_1`/`_2` × `.fastq[.gz]`/`.fq[.gz]`.

```bash
linezolid-amr folder -i input_dir/ -o results/
```

Produces `results/ALL_samples.summary_wide.csv` (one row per sample, Kleborate-style) and `results/ALL_samples.summary_long.csv` (one row per gene/mutation).

### All flags

```
linezolid-amr run / folder
  -a, --assembly       FASTA assembly                       (run only)
  -1, --r1             Forward reads                        (run only)
  -2, --r2             Reverse reads                        (run only)
  -i, --input          Folder with assemblies+reads         (folder only)
  -o, --outdir         Output directory                     (required)
  -s, --sample         Sample name                          (default: assembly stem)
  -O, --organism       Override MLST-inferred organism
  -t, --threads        CPUs (default: all available)
      --plus           AMRFinderPlus --plus (stress/virulence/biocide)
      --min-af 0.01    Min alt-allele frequency at 23S
      --min-depth 20   Min depth at 23S
      --skip-amrfinder
      --skip-rrna23s
```

---

## 📁 Output

```
results/sample/
├── amrfinder/amrfinder.tsv                    # all AMR/virulence hits
├── rrna23s/
│   ├── sample.23S.bam                          # sorted, indexed alignment
│   ├── sample.23S.vcf.gz                       # all 23S variants (bcftools)
│   └── sample.23S_lzd_pileup.tsv               # allele frequencies at LZD sites
├── sample.linezolid_amr.json                   # combined machine-readable report
├── sample.linezolid_amr.txt                    # plain-text summary
├── sample.summary_wide.csv                     # 1 row, AMR buckets as columns
└── sample.summary_long.csv                     # 1 row per gene/mutation
```

`folder` mode adds two cohort-level CSVs at the top of `outdir/`.

---

## 🧫 Worked example output

```
=== test ===
>> MLST / organism inference...
   organism: Enterococcus_faecium   MLST scheme: efaecium   ST: 17
>> Running AMRFinderPlus...
   17 hits, 3 linezolid-relevant
>> Running 23S rRNA analysis...
   11 positions; 1 with resistance allele

Linezolid resistance call: POSITIVE
```

23S pileup row:
```
ecoli_pos  ref  depth  counts        alt_alleles            is_resistance
2576       G    173    G=58; T=115   T:115:0.6647*          True   ← 4/6 operoni mutati
```

---

## 📚 References

All canonical 23S linezolid-resistance positions carry verified PMID citations in [`linezolid_amr/data/references/loci.json`](linezolid_amr/data/references/loci.json). See [`CITATION.cff`](CITATION.cff) for the full bibliography.

Key papers:
- Kloss et al. 1999 (foundational mutational mapping) — PMID 10556031
- Tsiodras et al. 2001 (first clinical G2576T) — PMID 11476839
- Long & Vester 2012 (comprehensive review) — PMID 22143525
- Long et al. 2006 (Cfr methyltransferase) — PMID 16801432
- Wang et al. 2015 (optrA discovery) — PMID 25977397
- Antonelli et al. 2018 (poxtA discovery) — PMID 29635422

---

## 📜 License

MIT — see [LICENSE](LICENSE).
