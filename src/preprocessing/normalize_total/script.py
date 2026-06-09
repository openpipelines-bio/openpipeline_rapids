import sys
import rapids_singlecell as rsc

## VIASH START
par = {
    "input": "work/d9/3adbd080e0de618d44b59b1ec81685/run.output.h5mu",
    "output": "output.h5mu",
    "target_sum": 10000,
    "modality": "rna",
    "exclude_highly_expressed": False,
    "max_fraction": 0.05,
}
meta = {"name": "normalize_total"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from anndata_io import read_modality, write_modality
from gpu import on_gpu

logger = setup_logger()

dat = read_modality(par, logger)

logger.info(par)

if par["input_layer"] and par["input_layer"] not in dat.layers.keys():
    raise ValueError(f"Input layer {par['input_layer']} not found in {par['modality']}")

if par["output_layer"]:
    source = dat.layers[par["input_layer"]] if par["input_layer"] else dat.X
    dat.layers[par["output_layer"]] = source.copy()

target_layer = par["output_layer"] or par["input_layer"]

with on_gpu(dat, logger, layer=target_layer):
    logger.info("Performing total normalization.")
    rsc.pp.normalize_total(
        dat,
        layer=target_layer,
        target_sum=par["target_sum"],
        exclude_highly_expressed=par["exclude_highly_expressed"],
        max_fraction=par["max_fraction"],
    )

write_modality(par, dat, logger)
