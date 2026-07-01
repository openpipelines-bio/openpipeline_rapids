#!/bin/bash

set -eo pipefail

REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

export NXF_SYNTAX_PARSER=v1

nextflow \
  run . \
  -main-script src/workflows/tools/leiden/test.nf \
  -entry test_wf \
  -profile docker,no_publish \
  -c src/workflows/utils/labels_ci.config \
  -c src/workflows/utils/workflow_tests.config
