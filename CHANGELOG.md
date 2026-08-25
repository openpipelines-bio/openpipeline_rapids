# openpipeline_rapids x.x.x

## BUG FIXES

* Replaced non-ASCII characters in configs and scripts with ASCII equivalents (PR #101).

# openpipeline_rapids v0.1.2

## BUG FIXES

* Passing an explicit empty string for `--layer`/`--input_layer` (a sentinel for "use .X")
  no longer crashes GPU preprocessing components with `KeyError: ''`; it is now normalized
  to `None` like an omitted argument (PR#100).

# openpipeline_rapids v0.1.1

## BUG FIXES

* `preprocessing/highly_variable_genes`: the boolean `highly_variable` column is now also
  included in the `.varm` output instead of being dropped while moving the per-gene HVG
  metrics out of `.var` (PR#98).

# openpipeline_rapids v0.1.0

## NEW FEATURES

Initial release of **OpenPipeline Rapids**: GPU-accelerated single-cell components built on
[rapids-singlecell](https://rapids-singlecell.readthedocs.io/), targeting x86_64 + CUDA 13.
Components wrap their rapids-singlecell equivalents and mirror the interfaces of their CPU
counterparts in [openpipeline](https://github.com/openpipelines-bio/openpipeline) so they can be
swapped in 1:1.

* Preprocessing components: `preprocessing/normalize_total`, `preprocessing/log1p`,
  `preprocessing/scale`, `preprocessing/regress_out`, `preprocessing/pca`,
  `preprocessing/neighbors`, `preprocessing/highly_variable_genes`,
  `preprocessing/calculate_qc_metrics`, `preprocessing/filter_cells`,
  `preprocessing/filter_genes`, `preprocessing/scrublet`, `preprocessing/harmony_integrate`,
  `preprocessing/bbknn`.

* Tools components: `tools/umap`, `tools/tsne`, `tools/leiden`.

* Squidpy components: `squidpy/spatial_autocorr`.

* Workflows: `workflows/rna/log_normalize` and `workflows/multiomics/neighbors_leiden_umap`,
  plus CPU/GPU wrapper workflows under `wrappers/` that switch between the in-repo GPU
  components and their upstream openpipeline CPU equivalents.
