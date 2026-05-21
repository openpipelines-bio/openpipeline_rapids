import sys
import rapids_singlecell as rsc
import mudata as mu

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
    "seed": None,
    "obsm_output": "X_pca",
    "varm_output": "pca_loadings",
    "uns_output": "pca_variance",
    "overwrite": False,
}
meta = {"name": "pca"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()

logger.info("Reading modality %s from %s", par["modality"], par["input"])
dat = mu.read_h5ad(par["input"], mod=par["modality"])
assert dat.var_names.is_unique, "The var_names of the input modality must be be unique."

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

# Verify output slots are free (or overwrite is allowed) before doing GPU work
check_exist = {
    "obsm_output": ("obsm", par["obsm_output"]),
    "varm_output": ("varm", par["varm_output"]),
    "uns_output": ("uns", par["uns_output"]),
}
for parameter_name, (field, key) in check_exist.items():
    if key in getattr(dat, field):
        if not par["overwrite"]:
            raise ValueError(
                f"Requested to create field {key} in .{field} for modality "
                f"{par['modality']}, but field already exists."
            )
        del getattr(dat, field)[key]

logger.info("Transferring data to GPU.")
rsc.get.anndata_to_GPU(dat, layer=par["layer"])

logger.info("Computing PCA.")
rsc.pp.pca(
    dat,
    n_comps=par["num_components"],
    layer=par["layer"],
    mask_var=mask_var,
    chunked=par["chunked"],
    chunk_size=par["chunk_size"],
    random_state=par["seed"],
)

logger.info("Transferring data back to CPU.")
rsc.get.anndata_to_CPU(dat)

# rapids-singlecell stores results under fixed keys ("X_pca", "PCs", "pca").
# Rename them to the requested output slots, mirroring the openpipeline CPU
# pca interface.
if par["obsm_output"] != "X_pca":
    dat.obsm[par["obsm_output"]] = dat.obsm.pop("X_pca")
if par["varm_output"] != "PCs":
    dat.varm[par["varm_output"]] = dat.varm.pop("PCs")
pca_uns = dat.uns.pop("pca")
dat.uns[par["uns_output"]] = {
    "variance": pca_uns["variance"],
    "variance_ratio": pca_uns["variance_ratio"],
}

logger.info(
    "Writing to file to %s with compression %s",
    par["output"],
    par["output_compression"],
)
write_h5ad_to_h5mu_with_compression(
    par["output"], par["input"], par["modality"], dat, par["output_compression"]
)
