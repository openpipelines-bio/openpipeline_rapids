#!/bin/bash

set -eo pipefail

# Regenerate the pbmc_1k_protein_v3 test resources consumed by this package.
#
# The RNA components and workflows in openpipeline_rapids read two files:
#   - pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu  (raw counts as MuData)
#   - pbmc_1k_protein_v3_mms.h5mu                          (multi-sample, dimred + integrated)
#
# Both are produced by driving the upstream openpipeline pipelines directly from
# GitHub, which ships the built Nextflow target (referencing the pinned prebuilt
# images) on its release tags, so there is no local checkout or build.

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the commands below are run from the root of the repository
cd "$REPO_ROOT"

# openpipeline release pinned in _viash.yaml (repositories: openpipeline)
OP_TAG=$(grep -A3 '^  - name: openpipeline$' _viash.yaml | grep 'tag:' | head -1 | awk '{print $2}')

ID=pbmc_1k_protein_v3
OUT=resources_test/$ID/$ID
DIR=$(dirname "$OUT")

[ -d "$DIR" ] || mkdir -p "$DIR"

# dataset page:
# https://www.10xgenomics.com/resources/datasets/1-k-pbm-cs-from-a-healthy-donor-gene-expression-and-cell-surface-protein-3-standard-3-0-0

# download metrics summary
wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_metrics_summary.csv \
  -O "${OUT}_metrics_summary.csv"

# download counts h5 file
wget https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5 \
  -O "${OUT}_filtered_feature_bc_matrix.h5"

# convert 10x h5 to h5mu (consumed file #1)
nextflow run openpipelines-bio/openpipeline \
  -r "$OP_TAG" \
  -main-script target/nextflow/convert/from_10xh5_to_h5mu/main.nf \
  -profile docker \
  -c src/workflows/utils/labels_ci.config \
  --id "$ID" \
  --input "${OUT}_filtered_feature_bc_matrix.h5" \
  --input_metrics_summary "${OUT}_metrics_summary.csv" \
  --output "$(basename "$OUT")_filtered_feature_bc_matrix.h5mu" \
  --publishDir "$DIR" \
  -resume

# run single sample
nextflow run openpipelines-bio/openpipeline \
  -r "$OP_TAG" \
  -main-script target/nextflow/workflows/rna/rna_singlesample/main.nf \
  -profile docker \
  -c src/workflows/utils/labels_ci.config \
  --id "${ID}_uss" \
  --input "${OUT}_filtered_feature_bc_matrix.h5mu" \
  --output "$(basename "$OUT")_uss.h5mu" \
  --publishDir "$DIR" \
  -resume

# add the sample ID to the mudata object
nextflow run openpipelines-bio/openpipeline \
  -r "$OP_TAG" \
  -main-script target/nextflow/metadata/add_id/main.nf \
  -profile docker \
  -c src/workflows/utils/labels_ci.config \
  --id "${ID}_uss" \
  --input "${OUT}_uss.h5mu" \
  --input_id "${ID}_uss" \
  --output "$(basename "$OUT")_uss_with_id.h5mu" \
  --output_compression "gzip" \
  --publishDir "$DIR" \
  -resume

# run multisample
nextflow run openpipelines-bio/openpipeline \
  -r "$OP_TAG" \
  -main-script target/nextflow/workflows/rna/rna_multisample/main.nf \
  -profile docker \
  -c src/workflows/utils/labels_ci.config \
  --id "${ID}_ums" \
  --input "${OUT}_uss_with_id.h5mu" \
  --output "$(basename "$OUT")_ums.h5mu" \
  --publishDir "$DIR" \
  -resume

# run dimensionality reduction
nextflow run openpipelines-bio/openpipeline \
  -r "$OP_TAG" \
  -main-script target/nextflow/workflows/multiomics/dimensionality_reduction/main.nf \
  -profile docker \
  -c src/workflows/utils/labels_ci.config \
  --id "${ID}_mms" \
  --input "${OUT}_ums.h5mu" \
  --output "$(basename "$OUT")_mms.h5mu" \
  --publishDir "$DIR" \
  --obs_covariates sample_id \
  -resume

# run integration (overwrites _mms.h5mu with the harmony + leiden result; consumed file #2)
nextflow run openpipelines-bio/openpipeline \
  -r "$OP_TAG" \
  -main-script target/nextflow/workflows/integration/harmony_leiden/main.nf \
  -profile docker \
  -c src/workflows/utils/labels_ci.config \
  --id "${ID}_mms_integration" \
  --input "${OUT}_mms.h5mu" \
  --output "$(basename "$OUT")_mms.h5mu" \
  --publishDir "$DIR" \
  --obs_covariates sample_id \
  -resume

# drop everything else (intermediates, raw downloads, Nextflow .state.yaml
# sidecars); keep only the two files consumed by the tests
BASE=$(basename "$OUT")
find "$DIR" -mindepth 1 \
  ! -name "${BASE}_filtered_feature_bc_matrix.h5mu" \
  ! -name "${BASE}_mms.h5mu" \
  -delete

# sync to this package's own bucket (drop --dryrun to perform the upload)
aws s3 sync \
  --profile di \
  "$DIR" \
  s3://openpipelines-bio/openpipeline_rapids/resources_test/$ID \
  --delete \
  --dryrun
