# openpipeline_rapids vx.x.x

## MINOR CHANGES

* Bumped `rapids-singlecell-cu13` from 0.15.0 to 0.18.1 (PR #103). Upstream changes that are
  visible in component output:
  * `preprocessing/filter_cells` and `preprocessing/filter_genes` only annotate the count column
    that matches the threshold used (`n_counts` for `--min_counts`/`--max_counts`, `n_genes` or
    `n_cells` otherwise), instead of both.
  * `preprocessing/harmony_integrate` now computes in float32, uses a k-means++ initialization
    and models multiple `--obs_covariates` separately, so seeded results differ from earlier
    versions.
  * `tools/umap` runs the correct number of epochs for inputs above 10,000 cells, which changes
    embeddings for large datasets.
  * `preprocessing/bbknn` connectivities now match the upstream bbknn package.
  * `squidpy/spatial_autocorr` normalizes Geary's C p-values like squidpy 1.8.2.

* `preprocessing/highly_variable_genes`: missing values in the boolean `highly_variable` and
  `highly_variable_intersection` columns are filled with `False` (PR #105).

## BUG FIXES

* `preprocessing/pca`: a `--var_input` column containing NA values is now rejected with a clear
  error. Previously `pd.NA` crashed the boolean cast and `NaN` was silently read as `True`,
  pulling those genes into the PCA (PR #105).

# openpipeline_rapids v0.1.3

## BUG FIXES

* Replaced non-ASCII characters in configs and scripts with ASCII equivalents (PR #101).

# openpipeline_rapids v0.1.2

## BUG FIXES

* Passing an explicit empty string for `--layer`/`--input_layer` (a sentinel for "use .X")
  no longer crashes GPU preprocessing components with `KeyError: ''`; it is now normalized
  to `None` like an omitted argument (PR #100).

# openpipeline_rapids v0.1.1

## BUG FIXES

* `preprocessing/highly_variable_genes`: the boolean `highly_variable` column is now also
  included in the `.varm` output instead of being dropped while moving the per-gene HVG
  metrics out of `.var` (PR #98).

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
