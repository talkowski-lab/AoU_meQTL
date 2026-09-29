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

Converts a VCF into a PLINK1 binary fileset (`.bed`/`.bim`/`.fam`), applying standard QC filters.

- Single task (`Plink2MakeBed`) runs `plink2 --vcf ... --maf ~{MinAF} --hwe ~{HWEPvalThreshold} --make-bed` inside the `bioinformatics` Docker image.
- `MinAF` (default `0.01`) drops variants below that minor allele frequency; `HWEPvalThreshold` (default `1e-6`) drops variants failing the Hardy-Weinberg exact test at that p-value.
- Docker image is built from **this repo** (`envs/Dockerfile.bioinformatics`, plink2 + bcftools/tabix) and published to Docker Hub as `<DOCKERHUB_USERNAME>/<repo-name>-bioinformatics`. The WDL selects the tag via the `ImageTag` input (defaults to `latest`; pass a 7-char commit SHA to pin a specific build).

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

`.github/workflows/bioinformatics-docker-image.yml` builds and pushes the bioinformatics image to Docker Hub on push/PR to `main`/`develop`, currently only when `envs/Dockerfile.bioinformatics` changes — add a `scripts/**` path filter once a task in this env starts using a script from `scripts/`. Images are named `<DOCKERHUB_USERNAME>/<repo-name>-bioinformatics` and tagged `latest` + the 7-char commit SHA. Auth uses the `DOCKERHUB_USERNAME` repo **variable** and the `DOCKERHUB_TOKEN` **secret**.

## Gotchas

- `plink2` has no Debian package, so the Dockerfile pulls the prebuilt binary from cog-genomics' `_latest` S3 alias. Pin `PLINK2_VERSION` (build arg) to a dated build instead for a reproducible image.
- `--maf`/`--hwe` in plink2 apply across the whole sample by default. If case/control-aware HWE filtering is ever needed (only filtering on controls), that requires `--hwe ... midp` semantics or splitting by phenotype first — not currently implemented.
- The Docker Hub image name (`ayenkin1871/aou-mqtl-analysis-bioinformatics` in `vcf_to_plink.wdl`) is hardcoded to match this repo's expected name/owner; update it if the repo is renamed or forked under a different Docker Hub account.
