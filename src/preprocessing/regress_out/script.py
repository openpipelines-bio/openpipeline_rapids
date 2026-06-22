import sys
import rapids_singlecell as rsc

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "input_layer": None,
    "output_layer": None,
    "obs_keys": ["total_counts"],
    "batchsize": None,
    "output_compression": None,
}
meta = {"name": "regress_out"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from anndata_io import read_modality, write_modality
from gpu import on_gpu

logger = setup_logger()

dat = read_modality(par["input"], par["modality"], logger)

logger.info(par)

if par["input_layer"] and par["input_layer"] not in dat.layers.keys():
    raise ValueError(f"Input layer {par['input_layer']} not found in {par['modality']}")

missing_keys = [key for key in par["obs_keys"] if key not in dat.obs.columns]
if missing_keys:
    raise ValueError(
        f"Columns {missing_keys} not found in .obs of modality {par['modality']}"
    )

if par["output_layer"]:
    source = dat.layers[par["input_layer"]] if par["input_layer"] else dat.X
    dat.layers[par["output_layer"]] = source.copy()

target_layer = par["output_layer"] or par["input_layer"]

with on_gpu(dat, logger, layer=target_layer):
    logger.info("Regressing out covariates.")
    rsc.pp.regress_out(
        dat,
        keys=par["obs_keys"],
        layer=target_layer,
        batchsize=par["batchsize"],
    )

write_modality(
    dat, par["output"], par["input"], par["modality"], par["output_compression"], logger
)
