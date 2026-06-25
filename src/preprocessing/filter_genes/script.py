import sys
import rapids_singlecell as rsc

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "min_counts": None,
    "min_cells": 3,
    "max_counts": None,
    "max_cells": None,
    "output_compression": None,
}
meta = {"name": "filter_genes"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from anndata_io import read_modality, write_modality
from gpu import on_gpu

logger = setup_logger()

dat = read_modality(par["input"], par["modality"], logger)

logger.info(par)

# rsc.pp.filter_genes accepts only a single threshold per call.
if (
    sum(
        par[threshold] is not None
        for threshold in ("min_counts", "min_cells", "max_counts", "max_cells")
    )
    != 1
):
    raise ValueError(
        "Exactly one of --min_counts, --min_cells, --max_counts, --max_cells "
        "must be set per call."
    )

with on_gpu(dat, logger):
    logger.info("Filtering genes.")
    rsc.pp.filter_genes(
        dat,
        min_counts=par["min_counts"],
        min_cells=par["min_cells"],
        max_counts=par["max_counts"],
        max_cells=par["max_cells"],
    )

write_modality(
    dat, par["output"], par["input"], par["modality"], par["output_compression"], logger
)
