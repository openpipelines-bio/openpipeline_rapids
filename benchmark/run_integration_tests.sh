#!/bin/bash

# One-shot GPU VM smoke test for the benchmark workflows.
#
# Verifies that both the GPU (process_gpu) and CPU (process_cpu) pipelines run
# end to end on the small xenium fixture before launching the real benchmark on
# Seqera. Run this on an x86_64 host with an NVIDIA GPU, Docker + the NVIDIA
# Container Toolkit, viash, nextflow and the AWS CLI installed.

set -eo pipefail

REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

# 1. Sync the xenium test fixture (public bucket, ~1.4 MB).
echo ">> Syncing xenium test data"
aws s3 sync \
  s3://openpipelines-bio/openpipeline_spatial/resources_test/xenium \
  resources_test/xenium \
  --no-sign-request

# 2. Build all components and workflows, building the rapids Docker images as
#    needed. The CPU OpenPipeline images are pulled from ghcr.io at run time.
echo ">> Building components and workflows"
viash ns build --parallel --setup cachedbuild

# 3. Run both integration tests. Each maps the `gpu` label to `--gpus all` via
#    labels_ci.config; the CPU pipeline simply ignores it.
echo ">> Running process_gpu integration test"
bash src/workflows/benchmark/process_gpu/integration_test.sh

echo ">> Running process_cpu integration test"
bash src/workflows/benchmark/process_cpu/integration_test.sh

echo ">> Both integration tests passed"
