#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

# releases pinned in _viash.yaml (repositories:)
SPATIAL_TAG=$(grep -A3 '^  - name: openpipeline_spatial$' _viash.yaml | grep 'tag:' | head -1 | awk '{print $2}')
OP_TAG=$(grep -A3 '^  - name: openpipeline$' _viash.yaml | grep 'tag:' | head -1 | awk '{print $2}')

SPATIAL_REPO=https://packages.viash-hub.com/vsh/openpipeline_spatial

DIR="$REPO_ROOT/resources_test/xenium"
ID="xenium_tiny"
OUT="$DIR/$ID"

# create tempdir for the raw download
MY_TEMP="${VIASH_TEMP:-/tmp}"
TMPDIR=$(mktemp -d "$MY_TEMP/$ID-XXXXXX")
function clean_up {
  [[ -d "$TMPDIR" ]] && rm -r "$TMPDIR"
}
trap clean_up EXIT

# shallow-clone the built openpipeline_spatial package (pinned tag) so its
# Nextflow modules can be run from the local checkout
SPATIAL_DIR="$TMPDIR/openpipeline_spatial"
git clone --depth 1 --branch "$SPATIAL_TAG" "$SPATIAL_REPO" "$SPATIAL_DIR"

# download the nf-core tiny Xenium dataset
if [ ! -d "$OUT" ]; then
    tiny_dataset="https://raw.githubusercontent.com/nf-core/test-datasets/spatialxe/Xenium_Prime_Mouse_Ileum_tiny_outs.tar.gz"
    wget "$tiny_dataset" -O "$TMPDIR/xenium_tiny.tar.gz"

    mkdir -p "$TMPDIR/xenium_tiny"
    tar -xzf "$TMPDIR/xenium_tiny.tar.gz" -C "$TMPDIR/xenium_tiny"
    mkdir -p "$OUT"
    mv "$TMPDIR/xenium_tiny/Xenium_Prime_Mouse_Ileum_tiny_outs/"* "$OUT/"
fi

# xenium -> spatialdata (.zarr)
rm -rf "$DIR/$ID.zarr"
nextflow run "$SPATIAL_DIR/target/nextflow/convert/from_xenium_to_spatialdata/main.nf" \
  -profile docker \
  -c src/workflows/utils/labels_ci.config \
  --id "$ID" \
  --input "$OUT" \
  --output "$ID.zarr" \
  --publishDir "$DIR" \
  -resume

# spatialdata -> h5mu
nextflow run "$SPATIAL_DIR/target/nextflow/convert/from_spatialdata_to_h5mu/main.nf" \
  -profile docker \
  -c src/workflows/utils/labels_ci.config \
  --id "$ID" \
  --input "$DIR/$ID.zarr" \
  --output "$ID.h5mu" \
  --publishDir "$DIR" \
  -resume

# spatial neighborhood graph
nextflow run "$SPATIAL_DIR/target/nextflow/neighbors/spatial_neighborhood_graph/main.nf" \
  -profile docker \
  -c src/workflows/utils/labels_ci.config \
  --id "$ID" \
  --input "$DIR/$ID.h5mu" \
  --output "$ID.neighbors.h5mu" \
  --publishDir "$DIR" \
  -resume

# QC (driven from the pinned openpipeline release)
cat > "$TMPDIR/qc.yaml" <<EOF
param_list:
  - id: $ID
    input: "$DIR/$ID.neighbors.h5mu"
var_name_mitochondrial_genes: mitochondrial
var_name_ribosomal_genes: ribosomal
output: '\$id.qc.h5mu'
output_compression: gzip
publish_dir: "$DIR"
EOF

nextflow run openpipelines-bio/openpipeline \
  -r "$OP_TAG" \
  -main-script target/nextflow/workflows/qc/qc/main.nf \
  -profile docker \
  -params-file "$TMPDIR/qc.yaml" \
  -c src/workflows/utils/labels_ci.config \
  -resume

# spatial neighborhood graph on the QC'd object
nextflow run "$SPATIAL_DIR/target/nextflow/neighbors/spatial_neighborhood_graph/main.nf" \
  -profile docker \
  -c src/workflows/utils/labels_ci.config \
  --id "$ID" \
  --input "$DIR/$ID.qc.h5mu" \
  --output "${ID}.qc.neighbors_unfiltered.h5mu" \
  --publishDir "$DIR" \
  -resume

# drop low-detection genes (consumed file): rsc.gr.spatial_autocorr's float32
# reduction yields nan/inf on constant or near-constant genes, which this tiny
# Xenium panel has several of. Mirrors the min_cells=3 filter already applied
# in src/squidpy/spatial_autocorr/test.py.
nextflow run . \
  -main-script target/nextflow/preprocessing/filter_genes/main.nf \
  -profile docker \
  -c src/workflows/utils/labels_ci.config \
  --id "$ID" \
  --input "$DIR/${ID}.qc.neighbors_unfiltered.h5mu" \
  --min_cells 3 \
  --output "${ID}.qc.neighbors.h5mu" \
  --publishDir "$DIR" \
  -resume

# drop everything else (raw download, intermediates, Nextflow .state.yaml
# sidecars); keep only the file consumed by the tests
find "$DIR" -mindepth 1 \
  ! -name "${ID}.qc.neighbors.h5mu" \
  -delete

# sync to this package's own bucket (drop --dryrun to perform the upload)
aws s3 sync \
  --profile di \
  "$DIR" \
  s3://openpipelines-bio/openpipeline_rapids/resources_test/xenium \
  --delete \
  --dryrun
