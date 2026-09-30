# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Purpose

Bioinformatics workflows for molecular QTL (mQTL) analysis in the All of Us (AoU) cohort. WDL workflows registered for Dockstore, run on Cromwell/Terra.

## Architecture

### Overall Organization
- workflows/ contains any WDLs that will be used for pipelines, subfolders and utility WDL files are allowed
- scripts/ contains any scripts that will be used within custom environments for pipelines, every dockerfile should copy the entire scripts/ folder
- envs/ contain Dockerfiles for custom environments used in workflows, named flat as `Dockerfile.<env-name>` (not one Dockerfile per subfolder)

All envs should have their own rule in .github/workflows so that they will be updated whenever the relevant dockerfile or any of scripts/* is changed.
All workflows should be registered in .dockstore.yml.

### 1. VCFToPlink (`workflows/vcf_to_plink.wdl`)

Converts a VCF into a PLINK2 fileset (`.pgen`/`.pvar`/`.psam`), applying standard QC filters.

- Single task (`Plink2MakePgen`) runs `plink2 --vcf ... --maf ~{MinAF} --hwe ~{HWEPvalThreshold} --make-pgen` inside the `bioinformatics` Docker image.
- `MinAF` (default `0.01`) drops variants below that minor allele frequency; `HWEPvalThreshold` (default `1e-6`) drops variants failing the Hardy-Weinberg exact test at that p-value.
- Docker image is built from **this repo** (`envs/Dockerfile.bioinformatics`) and published to Docker Hub as `<DOCKERHUB_USERNAME>/aou_meqtl-bioinformatics` (CI lowercases the repo name — `AoU_meQTL` — since Docker Hub image names must be lowercase). The WDL selects the tag via the `ImageTag` input (defaults to `latest`; pass a 7-char commit SHA to pin a specific build).

### 2. BinCpGMethylation (`workflows/bin_cpg_methylation.wdl` + `scripts/bin_cpg_methylation.py`)

Bins per-CpG methylation calls (a pb-CpG-tools bed) into fixed-size genomic windows and summarizes each window.

- Optional first task (`IntersectWithIntervals`): if `IntervalBed` or `IntervalString` is given, restricts `CpGBed` to those regions via `bedtools intersect -u` before binning (`IntervalBed` takes precedence if both are set). `IntervalString` is a 1-based inclusive region, e.g. `chr1:1000000-2000000` (samtools/tabix style), converted to 0-based BED internally.
- `BinCpGs` runs `bedtools makewindows -g ~{ChromSizes} -w ~{WindowSize}` to tile the genome, then `bedtools intersect -wa -wb` against the (possibly filtered) CpG bed, piping pairs into `scripts/bin_cpg_methylation.py` to aggregate per window: `num_cpgs`, `total_coverage`, `weighted_mean_methylation` (coverage-weighted), `unweighted_mean_methylation` (simple mean across CpGs). Windows with zero overlapping CpGs are dropped rather than emitted with NAs.
- `MethCol`/`CovCol` (defaults `4`/`6`) are 1-based column indices into `CpGBed`, assuming the standard pb-CpG-tools combined pileup bed (`chrom, start, end, modification_probability, haplotype, coverage`). Adjust these if the input bed has a different layout (e.g. "count" mode, which adds modified/unmodified count columns).
- `ChromSizes` is a standard 2-column `chrom<TAB>size` file (as produced by `cut -f1,2 ref.fa.fai` or UCSC `chrom.sizes`).
- Output bed has a `#`-prefixed header and is sorted by `chrom,start`.
- Docker image: same `envs/Dockerfile.bioinformatics` as VCFToPlink (adds `bedtools`/`python3`; the script is baked in via `COPY scripts/ /scripts/`).

## Common Commands

There is no build/test/lint tooling in this repo — it is WDL + a Dockerfile.

Build the bioinformatics Docker image locally (note: Docker context is repo root, Dockerfile path is explicit):
```
docker build -f envs/Dockerfile.bioinformatics -t aou-mqtl-analysis .
```

Validate / run WDLs locally (requires miniwdl or womtool/cromwell, not vendored):
```
miniwdl check workflows/vcf_to_plink.wdl
```

## CI

`.github/workflows/bioinformatics-docker-image.yml` builds and pushes the bioinformatics image to Docker Hub on push/PR to `main`/`develop`, when `scripts/**` or `envs/Dockerfile.bioinformatics` changes. Images are named `<DOCKERHUB_USERNAME>/aou_meqtl-bioinformatics` and tagged `latest` + the 7-char commit SHA. Auth uses the `DOCKERHUB_USERNAME` repo **variable** and the `DOCKERHUB_TOKEN` **secret**.

## Gotchas

- `plink2` has no Debian package, so the Dockerfile pulls the prebuilt binary from cog-genomics' `_latest` S3 alias. Pin `PLINK2_VERSION` (build arg) to a dated build instead for a reproducible image.
- `--maf`/`--hwe` in plink2 apply across the whole sample by default. If case/control-aware HWE filtering is ever needed (only filtering on controls), that requires `--hwe ... midp` semantics or splitting by phenotype first — not currently implemented.
- The repo name (`AoU_meQTL`) has uppercase letters, but Docker Hub image names must be lowercase — the CI workflow lowercases `REPO_NAME` before building the tag. The hardcoded Docker image string in `vcf_to_plink.wdl` (`ayenkin1871/aou_meqtl-bioinformatics`) must be kept in sync with whatever that lowercased name resolves to; update both if the repo is renamed or forked under a different Docker Hub account.
