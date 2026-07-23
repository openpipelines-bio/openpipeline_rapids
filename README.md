# OpenPipeline Rapids

GPU-accelerated single-cell analysis components for reproducible and large-scale single-cell processing using Viash and Nextflow.

OpenPipeline Rapids extends the [OpenPipeline](https://github.com/openpipelines-bio/openpipeline/) ecosystem with GPU-enabled components and workflows built on top of [rapids-singlecell](https://rapids-singlecell.readthedocs.io/). The structure of this package mirrors the rapids-singlecell Python package, and many components replicate their core OpenPipeline CPU equivalents with compatible interfaces so they can be swapped in 1:1 to accelerate existing pipelines on NVIDIA GPUs.

[![ViashHub](https://img.shields.io/badge/ViashHub-openpipeline_rapids-7a4baa.svg)](https://www.viash-hub.com/packages/openpipeline_rapids)
[![GitHub](https://img.shields.io/badge/GitHub-openpipelines--bio%2Fopenpipeline_rapids-blue.svg)](https://github.com/openpipelines-bio/openpipeline_rapids)
[![GitHub License](https://img.shields.io/github/license/openpipelines-bio/openpipeline_rapids.svg)](https://github.com/openpipelines-bio/openpipeline_rapids/blob/main/LICENSE)
[![GitHub Issues](https://img.shields.io/github/issues/openpipelines-bio/openpipeline_rapids.svg)](https://github.com/openpipelines-bio/openpipeline_rapids/issues)
[![Viash version](https://img.shields.io/badge/Viash-v0.9.7-blue.svg)](https://viash.io)

## Functionality

OpenPipeline Rapids provides GPU-accelerated implementations of common single-cell processing steps as standalone components with a standardized interface, along with wrappers that let you choose between the CPU (upstream OpenPipeline) and GPU (rapids-singlecell) implementation at runtime.

The following functionality is provided:

- [Preprocessing](https://www.viash-hub.com/packages/openpipeline_rapids/latest/components?search=preprocessing): QC metric calculation, cell and gene filtering, doublet detection (Scrublet), total-count normalization, log1p transformation, scaling, regressing out unwanted variation, highly variable gene selection, PCA, neighbor graph construction, and batch integration (Harmony, BBKNN).
- [Tools](https://www.viash-hub.com/packages/openpipeline_rapids/latest/components?search=tools): Clustering (Leiden) and dimensionality reduction / embedding (UMAP, t-SNE).
- [Spatial](https://www.viash-hub.com/packages/openpipeline_rapids/latest/components?search=squidpy): Spatial autocorrelation statistics (Squidpy) for spatial omics data.
- [Wrappers](https://www.viash-hub.com/packages/openpipeline_rapids/latest/components?search=wrappers): Drop-in CPU/GPU variants of the components and workflows above, exposing a `--device_type` switch so the same pipeline step can run on CPU or GPU without changing its interface.

## Extended functionality

This package only provides GPU-accelerated components; it is designed to work seamlessly with the core [OpenPipeline package](https://github.com/openpipelines-bio/openpipeline/). Because the components share compatible interfaces with their CPU counterparts, they can be dropped into existing OpenPipeline workflows to accelerate individual steps, and all core OpenPipeline workflows and components can be used in conjunction with the GPU-accelerated ones.

## Requirements

The components run on NVIDIA GPUs only. Execution requires:

- An NVIDIA GPU with the [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/) installed (the executable runner passes `--gpus all` to Docker).
- Components target x86_64 + CUDA 13. The base image is `nvidia/cuda:13.1.2-runtime-ubuntu24.04`, with the RAPIDS stack (`cudf`, `cuml`, `cugraph`, `cuvs`, `rapids-singlecell`) installed from NVIDIA's package index.

## Execution via CLI or Seqera Cloud

The openpipeline_rapids package is available via [Viash Hub](https://www.viash-hub.com/packages/openpipeline_rapids/latest/), where you can find instructions on how to run the workflows and individual components.

* Navigate to the [Viash Hub package page](https://www.viash-hub.com/packages/openpipeline_rapids/latest/), select the workflow or component you want to launch and click the `launch` button.
* Select the execution environment of choice (e.g. `Seqera Cloud`, `CLI` or `Executable`).
* Fill in the form with the required parameters and launch the workflow.

## Support

For issues specific to GPU-accelerated components, please use the [GitHub issues tracker](https://github.com/openpipelines-bio/openpipeline_rapids/issues). For general OpenPipeline questions, refer to the main [OpenPipeline documentation](https://openpipelines.bio/).
</content>
</invoke>
