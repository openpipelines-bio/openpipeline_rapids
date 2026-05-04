# OpenPipeline Rapids

[![License](https://img.shields.io/github/license/openpipelines-bio/openpipeline_rapids.svg)](https://github.com/openpipelines-bio/openpipeline_rapids/blob/main/LICENSE)

GPU-accelerated single-cell analysis components for the [OpenPipeline](https://github.com/openpipelines-bio/openpipeline/) ecosystem, built on [rapids-singlecell](https://rapids-singlecell.readthedocs.io/) and packaged with [Viash](https://viash.io) and [Nextflow](https://www.nextflow.io/).

This package mirrors selected components from the core OpenPipeline package and reimplements them on top of NVIDIA RAPIDS so they can run on CUDA-capable GPUs. The interfaces are kept compatible with their CPU counterparts so that GPU and CPU components can be swapped in the same workflow.
