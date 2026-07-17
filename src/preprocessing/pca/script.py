import sys
import numpy as np
import rapids_singlecell as rsc

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "layer": None,
    "var_input": None,
    "num_components": 25,
    "chunked": False,
    "chunk_size": None,
    "random_state": None,
    "obsm_output": "X_pca",
    "varm_output": "pca_loadings",
    "uns_output": "pca_variance",
    "overwrite": False,
}
meta = {"name": "pca"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from anndata_io import read_modality, write_modality
from gpu import on_gpu

logger = setup_logger()

dat = read_modality(par["input"], par["modality"], logger)

logger.info(par)

if par["layer"] and par["layer"] not in dat.layers.keys():
    raise ValueError(f"{par['layer']} was not found in modality {par['modality']}.")

if par["chunked"] and not par["chunk_size"]:
    raise ValueError(
        "Requested to perform an incremental PCA ('chunked'), but the chunk "
        "size is not set."
    )
if (
    par["chunked"]
    and par["num_components"]
    and par["chunk_size"] < par["num_components"]
):
    raise ValueError(
        f"The requested chunk size ({par['chunk_size']}) must not be smaller "
        f"than the number of components ({par['num_components']})"
    )

mask_var = None
if par["var_input"]:
    if par["var_input"] not in dat.var.columns:
        raise ValueError(
            f"Requested to use .var column {par['var_input']} as a selection "
            f"of genes to run the PCA on, but the column is not available "
            f"for modality {par['modality']}"
        )
    mask_var = par["var_input"]

# Fail fast if an output slot already exists and overwrite is not allowed.
# The slots themselves do not need clearing: the PCA call and the rename
# below assign into them unconditionally, overwriting any existing value.
check_exist = {
    "obsm_output": ("obsm", par["obsm_output"]),
    "varm_output": ("varm", par["varm_output"]),
    "uns_output": ("uns", par["uns_output"]),
}
for parameter_name, (field, key) in check_exist.items():
    if key in getattr(dat, field) and not par["overwrite"]:
        raise ValueError(
            f"Requested to create field {key} in .{field} for modality "
            f"{par['modality']}, but field already exists."
        )

# rapids-singlecell PCA cannot handle genes with zero total expression (they
# have zero variance) and fails with an opaque error deep in the call. Detect
# them up front and raise an actionable message instead.
matrix = dat.layers[par["layer"]] if par["layer"] else dat.X
selected_var_names = dat.var_names
if mask_var is not None:
    selected = dat.var[mask_var].to_numpy().astype(bool)
    matrix = matrix[:, selected]
    selected_var_names = selected_var_names[selected]
gene_totals = np.asarray(matrix.sum(axis=0)).ravel()
zero_genes = selected_var_names[gene_totals == 0].tolist()
if zero_genes:
    preview = ", ".join(zero_genes[:10])
    if len(zero_genes) > 10:
        preview += ", ..."
    raise ValueError(
        f"{len(zero_genes)} gene(s) selected for PCA have zero total expression "
        f"across all cells, which rapids-singlecell PCA cannot handle: "
        f"{preview}. Filter these genes out first (e.g. with `filter_genes`)."
    )

with on_gpu(dat, logger, layer=par["layer"]):
    logger.info("Computing PCA.")
    rsc.pp.pca(
        dat,
        n_comps=par["num_components"],
        layer=par["layer"],
        mask_var=mask_var,
        chunked=par["chunked"],
        chunk_size=par["chunk_size"],
        random_state=par["random_state"],
    )

# rapids-singlecell stores results under fixed keys ("X_pca", "PCs", "pca").
# We rename them rather than using the `key_added` argument because that sets
# a single shared key for all three slots, whereas this component exposes
# independent obsm/varm/uns slot names.
if par["obsm_output"] != "X_pca":
    dat.obsm[par["obsm_output"]] = dat.obsm.pop("X_pca")
if par["varm_output"] != "PCs":
    dat.varm[par["varm_output"]] = dat.varm.pop("PCs")
pca_uns = dat.uns.pop("pca")
dat.uns[par["uns_output"]] = {
    "variance": pca_uns["variance"],
    "variance_ratio": pca_uns["variance_ratio"],
}

write_modality(
    dat, par["output"], par["input"], par["modality"], par["output_compression"], logger
)
