import sys
import rapids_singlecell as rsc

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "layer": None,
    "batch_key": None,
    "sim_doublet_ratio": 2.0,
    "expected_doublet_rate": 0.05,
    "stdev_doublet_rate": 0.02,
    "n_prin_comps": 30,
    "n_neighbors": None,
    "threshold": None,
    "random_state": 0,
    "output_compression": None,
}
meta = {"name": "scrublet"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from anndata_io import read_modality, write_modality
from gpu import on_gpu

logger = setup_logger()

dat = read_modality(par["input"], par["modality"], logger)

logger.info(par)

if par["layer"] and par["layer"] not in dat.layers.keys():
    raise ValueError(f"Layer {par['layer']} not found in modality {par['modality']}")

if par["batch_key"] and par["batch_key"] not in dat.obs.columns:
    raise ValueError(
        f"Batch key {par['batch_key']} not found in .obs of modality {par['modality']}"
    )

# rsc.pp.scrublet always reads from .X. When --layer is given, stage that layer
# into .X for the call and restore the original .X afterward.
original_x = None
if par["layer"]:
    original_x = dat.X
    dat.X = dat.layers[par["layer"]]

with on_gpu(dat, logger):
    logger.info("Predicting doublets with Scrublet.")
    rsc.pp.scrublet(
        dat,
        batch_key=par["batch_key"],
        sim_doublet_ratio=par["sim_doublet_ratio"],
        expected_doublet_rate=par["expected_doublet_rate"],
        stdev_doublet_rate=par["stdev_doublet_rate"],
        n_prin_comps=par["n_prin_comps"],
        n_neighbors=par["n_neighbors"],
        threshold=par["threshold"],
        random_state=par["random_state"],
    )

if par["layer"]:
    dat.X = original_x

write_modality(
    dat, par["output"], par["input"], par["modality"], par["output_compression"], logger
)
