import sys
import rapids_singlecell as rsc

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "obsm_input": "X_pca",
    "uns_output": "neighbors",
    "obsp_distances": "distances",
    "obsp_connectivities": "connectivities",
    "num_neighbors": 15,
    "n_pcs": None,
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
from anndata_io import read_modality, write_modality

logger = setup_logger()

dat = read_modality(par["input"], par["modality"], logger)

logger.info(par)

if par["obsm_input"] not in dat.obsm.keys():
    raise ValueError(
        f"obsm slot '{par['obsm_input']}' not found in modality '{par['modality']}'."
    )

# rsc.pp.neighbors writes to .uns and .obsp (no matrix transform), and reads
# the embedding directly from .obsm[use_rep]. No anndata_to_GPU pre-stage is
# required — the function transfers the obsm representation to GPU internally.

logger.info("Computing a neighborhood graph.")
rsc.pp.neighbors(
    dat,
    n_neighbors=par["num_neighbors"],
    n_pcs=par["n_pcs"],
    use_rep=par["obsm_input"],
    random_state=par["seed"],
    algorithm=par["algorithm"],
    metric=par["metric"],
    method=par["method"],
)

# rsc.pp.neighbors writes to the default slots (.uns["neighbors"],
# .obsp["distances"], .obsp["connectivities"]). Move them to the requested
# slots so the output matches the configured argument names.
neighbors_uns = dat.uns.pop("neighbors")
neighbors_uns["distances_key"] = par["obsp_distances"]
neighbors_uns["connectivities_key"] = par["obsp_connectivities"]
dat.uns[par["uns_output"]] = neighbors_uns
dat.obsp[par["obsp_distances"]] = dat.obsp.pop("distances")
dat.obsp[par["obsp_connectivities"]] = dat.obsp.pop("connectivities")

write_modality(
    dat, par["output"], par["input"], par["modality"], par["output_compression"], logger
)
