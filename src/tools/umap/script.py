import sys
import rapids_singlecell as rsc

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "uns_neighbors": "neighbors",
    "obsm_output": "X_umap",
    "min_dist": 0.5,
    "spread": 1.0,
    "num_components": 2,
    "max_iter": None,
    "alpha": 1.0,
    "negative_sample_rate": 5,
    "init_pos": "spectral",
    "random_state": 0,
}
meta = {"name": "umap"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from anndata_io import read_modality, write_modality

logger = setup_logger()

dat = read_modality(par["input"], par["modality"], logger)

logger.info(par)

if par["uns_neighbors"] not in dat.uns:
    raise ValueError(
        f"'{par['uns_neighbors']}' was not found in .mod['{par['modality']}'].uns."
    )

logger.info("Computing UMAP for modality '%s'.", par["modality"])
rsc.tl.umap(
    dat,
    min_dist=par["min_dist"],
    spread=par["spread"],
    n_components=par["num_components"],
    maxiter=par["max_iter"],
    alpha=par["alpha"],
    negative_sample_rate=par["negative_sample_rate"],
    init_pos=par["init_pos"],
    random_state=par["random_state"],
    key_added=par["obsm_output"],
    neighbors_key=par["uns_neighbors"],
)

write_modality(
    dat, par["output"], par["input"], par["modality"], par["output_compression"], logger
)
