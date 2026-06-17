# GPU vs CPU bulk benchmark

Two Viash/Nextflow workflows chain every analytical step in this repository into
a single run so the full GPU rapids-singlecell stack can be timed against the
equivalent CPU OpenPipeline stack on the same input.

| Step | `process_gpu` (this repo, GPU) | `process_cpu` (OpenPipeline, CPU) |
|------|--------------------------------|-----------------------------------|
| normalize + log1p | `workflows/rna/log_normalize` | `workflows/rna/log_normalize` |
| highly variable genes | `preprocessing/highly_variable_genes` | `feature_annotation/highly_variable_features_scanpy` |
| PCA | `preprocessing/pca` | `dimred/pca` |
| neighbors + leiden + umap | `workflows/multiomics/neighbors_leiden_umap` | `workflows/multiomics/neighbors_leiden_umap` |
| spatial autocorrelation | `squidpy/spatial_autocorr` | `feature_annotation/spatial_autocorr` |

These workflows are for timing only; they are not meant to produce meaningful
biological results.

## Input requirements

A spatial `.h5mu` with:
- raw counts in `.X` of the processed modality (default `rna`), and
- a spatial neighbourhood graph in `.obsp["spatial_connectivities"]` (consumed
  by the `spatial_autocorr` step).

## Smoke test on a GPU VM first

Before spending time on Seqera, confirm both pipelines run end to end on the
small xenium fixture (23 cells, already carrying `obsp["spatial_connectivities"]`).
On an x86_64 host with an NVIDIA GPU, Docker + the NVIDIA Container Toolkit,
`viash`, `nextflow` and the AWS CLI:

```bash
bash benchmark/run_integration_tests.sh
```

This syncs the test data, builds the rapids Docker images
(`viash ns build --setup cachedbuild`), and runs both
`integration_test.sh` scripts. The `gpu` label is mapped to `--gpus all` via
`labels_ci.config`; the CPU pipeline ignores it. The OpenPipeline CPU images are
pulled from ghcr.io at run time.

Note this is a correctness smoke test only -- `labels_ci.config` pins resources
and `maxForks = 1`, so the timings here are not representative. Use Seqera with
real resources and a real dataset for the actual benchmark.

## Running on Seqera

Point the pipeline at this repository, set the main script to the built
workflow, and supply the matching params file:

```bash
# GPU run (requires a GPU-mapped compute environment, see below)
nextflow run <this-repo> \
  -main-script target/nextflow/workflows/benchmark/process_gpu/main.nf \
  -profile docker \
  -params-file benchmark/process_gpu.params.yaml

# CPU run (CPU-only compute environment)
nextflow run <this-repo> \
  -main-script target/nextflow/workflows/benchmark/process_cpu/main.nf \
  -profile docker \
  -params-file benchmark/process_cpu.params.yaml
```

Fill in the `s3://YOUR-BUCKET/...` placeholders in each params file first, and
use the **same input** for both runs.

## GPU label mapping

Every component in `process_gpu` carries the inert `gpu` label. The compute
environment must map it to a real GPU directive or the run fails with no CUDA
device visible:

```groovy
withLabel: gpu { accelerator = 1 }                    // AWS/Google Batch, Kubernetes
withLabel: gpu { containerOptions = '--gpus all' }    // single-node Docker
```

See [`../src/workflows/utils/labels.config`](../src/workflows/utils/labels.config)
for the full set of resource labels.

## Comparing runtimes

Both workflows enable Nextflow's execution `trace` (via the resource configs).
Compare wall-clock per process and total across the two runs. For a fair
comparison, equalise CPU counts between the two compute environments, since the
CPU components scale with available cores while the GPU components do not.
