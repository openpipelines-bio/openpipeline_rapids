import sys
import rapids_singlecell as rsc

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "input_layer": None,
    "output_layer": None,
    "zero_center": True,
    "max_value": None,
    "output_compression": None,
}
meta = {"name": "scale"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from anndata_io import read_modality, write_modality
from gpu import on_gpu

logger = setup_logger()

dat = read_modality(par["input"], par["modality"], logger)

logger.info(par)

# An empty string is a sentinel for "use .X"; normalize to None.
if not par["input_layer"]:
    par["input_layer"] = None

if par["input_layer"] and par["input_layer"] not in dat.layers.keys():
    raise ValueError(f"Input layer {par['input_layer']} not found in {par['modality']}")

# Pre-stage the source matrix into the output layer so rsc.pp.scale can scale
# it in place there, leaving the input (.X or --input_layer) untouched.
if par["output_layer"]:
    source = dat.layers[par["input_layer"]] if par["input_layer"] else dat.X
    dat.layers[par["output_layer"]] = source.copy()

target_layer = par["output_layer"] or par["input_layer"]

with on_gpu(dat, logger, layer=target_layer):
    logger.info("Scaling data.")
    rsc.pp.scale(
        dat,
        layer=target_layer,
        zero_center=par["zero_center"],
        max_value=par["max_value"],
    )

write_modality(
    dat, par["output"], par["input"], par["modality"], par["output_compression"], logger
)
