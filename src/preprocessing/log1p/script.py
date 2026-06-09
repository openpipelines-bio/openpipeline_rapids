import sys
import rapids_singlecell as rsc

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "input_layer": None,
    "input_obsm": None,
    "output_layer": None,
    "output_obsm": None,
    "base": None,
}
meta = {"name": "log1p"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from anndata_io import read_modality, write_modality
from gpu import on_gpu

logger = setup_logger()

dat = read_modality(par, logger)

logger.info(par)

if par["input_layer"] and par["input_obsm"]:
    raise ValueError("--input_layer and --input_obsm are mutually exclusive.")
if par["output_layer"] and par["output_obsm"]:
    raise ValueError("--output_layer and --output_obsm are mutually exclusive.")
if par["output_obsm"] and not par["input_obsm"]:
    raise ValueError("--output_obsm requires --input_obsm to be set.")

if par["input_layer"] and par["input_layer"] not in dat.layers.keys():
    raise ValueError(f"Input layer {par['input_layer']} not found in {par['modality']}")
if par["input_obsm"] and par["input_obsm"] not in dat.obsm.keys():
    raise ValueError(f"Input obsm {par['input_obsm']} not found in {par['modality']}")

if par["output_layer"]:
    source = dat.layers[par["input_layer"]] if par["input_layer"] else dat.X
    dat.layers[par["output_layer"]] = source.copy()

if par["output_obsm"]:
    dat.obsm[par["output_obsm"]] = dat.obsm[par["input_obsm"]].copy()

target_layer = par["output_layer"] or par["input_layer"]
target_obsm = par["output_obsm"] or par["input_obsm"]

# rapids-singlecell modifies the target slot in place. Warn loudly when this
# means overwriting a named input that the user didn't explicitly copy.
if par["input_layer"] and not par["output_layer"]:
    logger.warning(
        "No --output_layer set; the transformation will be applied in-place "
        "to .layers['%s']. Set --output_layer to preserve the input layer.",
        par["input_layer"],
    )
if par["input_obsm"] and not par["output_obsm"]:
    logger.warning(
        "No --output_obsm set; the transformation will be applied in-place "
        "to .obsm['%s']. Set --output_obsm to preserve the input entry.",
        par["input_obsm"],
    )

# anndata_to_GPU has no obsm-specific argument; convert_all moves .obsm too.
gpu_kwargs = {"convert_all": True} if target_obsm else {"layer": target_layer}
with on_gpu(dat, logger, **gpu_kwargs):
    logger.info("Performing log1p transformation.")
    rsc.pp.log1p(
        dat,
        base=par["base"],
        layer=target_layer if not target_obsm else None,
        obsm=target_obsm,
    )

write_modality(par, dat, logger)
