import sys
import rapids_singlecell as rsc

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "layer": None,
    "qc_vars": None,
    "expr_type": "counts",
    "var_type": "genes",
    "log1p": True,
    "output_compression": None,
}
meta = {"name": "calculate_qc_metrics"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from anndata_io import read_modality, write_modality
from gpu import on_gpu

logger = setup_logger()

dat = read_modality(par["input"], par["modality"], logger)

logger.info(par)

if par["layer"] and par["layer"] not in dat.layers.keys():
    raise ValueError(
        f"Layer '{par['layer']}' not found in modality '{par['modality']}'."
    )

qc_vars = par["qc_vars"]
if qc_vars is not None and len(qc_vars) == 0:
    qc_vars = None

if qc_vars is not None:
    missing = [name for name in qc_vars if name not in dat.var.columns]
    if missing:
        raise ValueError(
            f"qc_vars {missing} not found in .var of modality '{par['modality']}'."
        )

with on_gpu(dat, logger, layer=par["layer"]):
    logger.info("Calculating QC metrics.")
    rsc.pp.calculate_qc_metrics(
        dat,
        expr_type=par["expr_type"],
        var_type=par["var_type"],
        qc_vars=qc_vars,
        log1p=par["log1p"],
        layer=par["layer"],
    )

write_modality(
    dat, par["output"], par["input"], par["modality"], par["output_compression"], logger
)
