import sys
import rapids_singlecell as rsc
import mudata as mu

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "obsm_input": "X_pca",
    "uns_output": "neighbors",
    "num_neighbors": 15,
    "metric": "euclidean",
    "algorithm": "brute",
    "method": "umap",
    "seed": 0,
    "output_compression": None,
}
meta = {"name": "neighbors"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()

logger.info("Reading modality %s from %s", par["modality"], par["input"])
dat = mu.read_h5ad(par["input"], mod=par["modality"])
assert dat.var_names.is_unique, "The var_names of the input modality must be be unique."

logger.info(par)

if par["obsm_input"] not in dat.obsm.keys():
    raise ValueError(
        f"obsm slot '{par['obsm_input']}' not found in modality '{par['modality']}'."
    )

# rsc.pp.neighbors writes to .uns and .obsp (no matrix transform), and reads
# the embedding directly from .obsm[use_rep]. No anndata_to_GPU pre-stage is
# required — the function transfers the obsm representation to GPU internally.

# rsc.pp.neighbors uses None to mean "default key", so we only pass key_added
# when the user requested a non-default uns slot.
key_added = par["uns_output"] if par["uns_output"] != "neighbors" else None

logger.info("Computing a neighborhood graph.")
rsc.pp.neighbors(
    dat,
    n_neighbors=par["num_neighbors"],
    use_rep=par["obsm_input"],
    random_state=par["seed"],
    algorithm=par["algorithm"],
    metric=par["metric"],
    method=par["method"],
    key_added=key_added,
)

logger.info(
    "Writing to file to %s with compression %s",
    par["output"],
    par["output_compression"],
)
write_h5ad_to_h5mu_with_compression(
    par["output"], par["input"], par["modality"], dat, par["output_compression"]
)
