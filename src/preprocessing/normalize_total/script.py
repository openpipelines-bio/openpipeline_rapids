import sys
import rapids_singlecell as rsc
import mudata as mu

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
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()

logger.info("Reading modality %s from %s", par["modality"], par["input"])
dat = mu.read_h5ad(par["input"], mod=par["modality"])
assert dat.var_names.is_unique, "The var_names of the input modality must be be unique."

logger.info(par)

if par["input_layer"] and par["input_layer"] not in dat.layers.keys():
    raise ValueError(f"Input layer {par['input_layer']} not found in {par['modality']}")

if par["output_layer"]:
    source = dat.layers[par["input_layer"]] if par["input_layer"] else dat.X
    dat.layers[par["output_layer"]] = source.copy()

target_layer = par["output_layer"] or par["input_layer"]

logger.info("Transferring data to GPU.")
rsc.get.anndata_to_GPU(dat, layer=target_layer)

logger.info("Performing total normalization.")
rsc.pp.normalize_total(
    dat,
    layer=target_layer,
    target_sum=par["target_sum"],
    exclude_highly_expressed=par["exclude_highly_expressed"],
    max_fraction=par["max_fraction"],
)

logger.info("Transferring data back to CPU.")
rsc.get.anndata_to_CPU(dat)

logger.info(
    "Writing to file to %s with compression %s",
    par["output"],
    par["output_compression"],
)
write_h5ad_to_h5mu_with_compression(
    par["output"], par["input"], par["modality"], dat, par["output_compression"]
)
