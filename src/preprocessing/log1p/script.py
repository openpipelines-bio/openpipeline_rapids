import sys
import rapids_singlecell as rsc
import mudata as mu

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "base": None,
}
meta = {"name": "log1p"}
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

logger.info("Transferring data to GPU.")
rsc.get.anndata_to_GPU(dat, layer=par["input_layer"])

logger.info("Performing log1p transformation.")
output_data = rsc.pp.log1p(
    dat,
    base=par["base"],
    layer=par["input_layer"],
    copy=True if par["output_layer"] else False,
)

logger.info("Transferring data back to CPU.")
rsc.get.anndata_to_CPU(dat)
if output_data:
    # separate copy returned by rsc.pp.log1p when copy=True
    rsc.get.anndata_to_CPU(output_data)
    result = (
        output_data.X
        if not par["input_layer"]
        else output_data.layers[par["input_layer"]]
    )
    dat.layers[par["output_layer"]] = result

logger.info(
    "Writing to file to %s with compression %s",
    par["output"],
    par["output_compression"],
)
write_h5ad_to_h5mu_with_compression(
    par["output"], par["input"], par["modality"], dat, par["output_compression"]
)
