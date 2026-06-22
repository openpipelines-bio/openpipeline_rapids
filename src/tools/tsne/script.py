import sys
import rapids_singlecell as rsc

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "obsm_input": "X_pca",
    "n_pcs": None,
    "obsm_output": "X_tsne",
    "perplexity": 30,
    "early_exaggeration": 12,
    "learning_rate": 200,
    "metric": "euclidean",
    "method": "barnes_hut",
    "output_compression": None,
}
meta = {"name": "tsne"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from anndata_io import read_modality, write_modality

logger = setup_logger()

dat = read_modality(par["input"], par["modality"], logger)

logger.info(par)

if par["obsm_input"] not in dat.obsm:
    raise ValueError(
        f"'{par['obsm_input']}' was not found in .mod['{par['modality']}'].obsm."
    )

logger.info("Computing t-SNE for modality '%s'.", par["modality"])
rsc.tl.tsne(
    dat,
    n_pcs=par["n_pcs"],
    use_rep=par["obsm_input"],
    perplexity=par["perplexity"],
    early_exaggeration=par["early_exaggeration"],
    learning_rate=par["learning_rate"],
    method=par["method"],
    metric=par["metric"],
    key_added=par["obsm_output"],
)

write_modality(
    dat, par["output"], par["input"], par["modality"], par["output_compression"], logger
)
