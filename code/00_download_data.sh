#!/usr/bin/env bash
# M0: Real data download (reproducible record).
# Network policy in this environment is an allowlist: UCSC Xena S3 hubs are
# reachable; CGGA (cgga.org.cn) and NCBI/GEO are NOT (HTTP 403).
# TCGA-GBM discovery cohort is therefore taken from UCSC Xena (TCGA Hub, legacy
# RNA-seqV2 = log2(norm_count+1)), which is the reachable equivalent of the
# GDC HTSeq matrix requested.
set -euo pipefail
cd "$(dirname "$0")/../data/raw"

XENA="https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download"

# Expression: IlluminaHiSeq RNAseqV2, gene-level, log2(norm_count+1)
curl -fL --retry 4 --retry-delay 2 -o TCGA_GBM_HiSeqV2.gz \
  "${XENA}/TCGA.GBM.sampleMap%2FHiSeqV2.gz"

# Curated survival (Liu et al. 2018 pan-cancer)
curl -fL --retry 4 --retry-delay 2 -o GBM_survival.txt \
  "${XENA}/survival%2FGBM_survival.txt"

# Clinical / phenotype matrix
curl -fL --retry 4 --retry-delay 2 -o TCGA_GBM_clinicalMatrix \
  "${XENA}/TCGA.GBM.sampleMap%2FGBM_clinicalMatrix"

echo "Downloaded:"
ls -la TCGA_GBM_HiSeqV2.gz GBM_survival.txt TCGA_GBM_clinicalMatrix
