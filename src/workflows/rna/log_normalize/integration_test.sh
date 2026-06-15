#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# pin Nextflow to a version that uses the legacy config syntax parser
export NXF_VER=25.10.6

nextflow \
  run . \
  -main-script src/workflows/rna/log_normalize/test.nf \
  -entry test_wf \
  -profile docker,no_publish \
  -c src/workflows/utils/labels_ci.config \
  -c src/workflows/utils/workflow_tests.config
