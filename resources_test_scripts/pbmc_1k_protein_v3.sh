#!/bin/bash

set -eo pipefail

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

# process the sample into a multi-sample object with dimensionality reduction;
# process_samples runs the full singlesample -> multisample -> dimred pipeline,
# producing the X_pca the components read
nextflow run openpipelines-bio/openpipeline \
  -r "$OP_TAG" \
  -main-script target/nextflow/workflows/multiomics/process_samples/main.nf \
  -profile docker \
  -c src/workflows/utils/labels_ci.config \
  --id "${ID}_mms" \
  --input "${OUT}_filtered_feature_bc_matrix.h5mu" \
  --output "$(basename "$OUT")_processed.h5mu" \
  --publishDir "$DIR" \
  -resume

# integrate (harmony) + cluster (leiden) -> consumed file #2; adds the
# harmony_integration_leiden_1.0 obs column and X_pca_integrated obsm that the
# integration wrapper workflow tests consume
nextflow run openpipelines-bio/openpipeline \
  -r "$OP_TAG" \
  -main-script target/nextflow/workflows/integration/harmony_leiden/main.nf \
  -profile docker \
  -c src/workflows/utils/labels_ci.config \
  --id "${ID}_mms_integration" \
  --input "${OUT}_processed.h5mu" \
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
