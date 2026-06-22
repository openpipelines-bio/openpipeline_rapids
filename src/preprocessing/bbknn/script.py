import sys
import rapids_singlecell as rsc

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "obsm_input": "X_pca",
    "batch_key": "batch",
    "uns_output": "neighbors",
    "obsp_distances": "distances",
    "obsp_connectivities": "connectivities",
    "neighbors_within_batch": 3,
    "n_pcs": None,
    "metric": "euclidean",
    "algorithm": "brute",
    "trim": None,
    "random_state": 0,
    "output_compression": None,
}
meta = {"name": "bbknn"}
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

if par["batch_key"] not in dat.obs.columns:
    raise ValueError(
        f"batch_key '{par['batch_key']}' not found in .obs of modality '{par['modality']}'."
    )

# rsc.pp.bbknn writes to .uns and .obsp (no matrix transform), and reads the
# embedding directly from .obsm[use_rep]. No anndata_to_GPU pre-stage is
# required - the function transfers the obsm representation to GPU internally.

logger.info("Computing a batch-balanced neighborhood graph.")
rsc.pp.bbknn(
    dat,
    neighbors_within_batch=par["neighbors_within_batch"],
    n_pcs=par["n_pcs"],
    batch_key=par["batch_key"],
    use_rep=par["obsm_input"],
    random_state=par["random_state"],
    algorithm=par["algorithm"],
    metric=par["metric"],
    trim=par["trim"],
)

# rsc.pp.bbknn writes to the default slots (.uns["neighbors"],
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
