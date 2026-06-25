import sys
import rapids_singlecell as rsc

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "layer": None,
    "min_counts": None,
    "min_genes": 50,
    "max_counts": None,
    "max_genes": None,
    "output_compression": None,
}
meta = {"name": "filter_cells"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from anndata_io import read_modality, write_modality
from gpu import on_gpu

logger = setup_logger()

dat = read_modality(par["input"], par["modality"], logger)

logger.info(par)

# rsc.pp.filter_cells accepts only a single threshold per call.
if (
    sum(
        par[threshold] is not None
        for threshold in ("min_counts", "min_genes", "max_counts", "max_genes")
    )
    != 1
):
    raise ValueError(
        "Exactly one of --min_counts, --min_genes, --max_counts, --max_genes "
        "must be set per call."
    )

if par["layer"] and par["layer"] not in dat.layers.keys():
    raise ValueError(f"Layer {par['layer']} not found in modality {par['modality']}")

# rsc.pp.filter_cells filters on .X. When --layer is given, swap it into .X for
# the call; filtering subsets .X and all layers consistently, so swapping back
# afterward leaves .X as the (subsetted) original matrix and the layer intact.
if par["layer"]:
    dat.X, dat.layers[par["layer"]] = dat.layers[par["layer"]], dat.X

# convert_all casts every matrix (including the swapped-in layer) to a
# GPU-compatible float dtype; without --layer only .X needs moving.
gpu_kwargs = {"convert_all": True} if par["layer"] else {}
with on_gpu(dat, logger, **gpu_kwargs):
    logger.info("Filtering cells.")
    rsc.pp.filter_cells(
        dat,
        min_counts=par["min_counts"],
        min_genes=par["min_genes"],
        max_counts=par["max_counts"],
        max_genes=par["max_genes"],
    )

if par["layer"]:
    dat.X, dat.layers[par["layer"]] = dat.layers[par["layer"]], dat.X

write_modality(
    dat, par["output"], par["input"], par["modality"], par["output_compression"], logger
)
