# OpenPipeline Rapids

GPU-accelerated single-cell analysis components for the [OpenPipeline](https://github.com/openpipelines-bio/openpipeline/) ecosystem, built on [rapids-singlecell](https://rapids-singlecell.readthedocs.io/) and packaged with [Viash](https://viash.io) and [Nextflow](https://www.nextflow.io/).

The structure of this package mirrors the rapids-singlecell Python package with components implementing rapids-singlecell functions.
Some core OpenPipeline components are replicated using small workflows with compatible interfaces allowing them to be swapped for the equivalent CPU component.

## Running on a GPU compute environment

Every component carries the Nextflow `gpu` label, but the label is intentionally inert by default because the directive that exposes a GPU depends on the executor. Your compute environment (e.g. on Seqera Platform) must define the mapping, otherwise the label is silently ignored and components fail at runtime with no CUDA device visible. Add one of the following to the compute environment's Nextflow config:

```groovy
// Local / single-node Docker executor:
withLabel: gpu { containerOptions = '--gpus all' }

// AWS/Google Batch or Kubernetes executors (requires a GPU-enabled queue):
withLabel: gpu { accelerator = 1 }
```

See [`src/workflows/utils/labels.config`](src/workflows/utils/labels.config) for the full set of resource labels.
