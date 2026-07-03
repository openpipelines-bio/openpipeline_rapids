# OpenPipeline Rapids

GPU-accelerated single-cell analysis components for the [OpenPipeline](https://github.com/openpipelines-bio/openpipeline/) ecosystem, built on [rapids-singlecell](https://rapids-singlecell.readthedocs.io/) and packaged with [Viash](https://viash.io) and [Nextflow](https://www.nextflow.io/).

The structure of this package mirrors the rapids-singlecell Python package with components implementing rapids-singlecell functions.
Some core OpenPipeline components are replicated using small workflows with compatible interfaces allowing them to be swapped for the equivalent CPU component.
