#!/bin/bash

set -eo pipefail

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

export NXF_SYNTAX_PARSER=v1

nextflow \
  run . \
  -main-script src/wrappers/squidpy/spatial_autocorr/test.nf \
  -entry test_wf \
  -profile docker,no_publish \
  -c src/workflows/utils/labels_ci.config \
  -c src/workflows/utils/workflow_tests.config
