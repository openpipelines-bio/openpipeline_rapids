import sys
import rapids_singlecell as rsc

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "obsm_input": "X_pca",
    "obs_covariates": ["batch"],
    "obsm_output": "X_pca_harmony",
    "theta": 2.0,
    "n_clusters": None,
    "max_iter_harmony": 10,
    "random_state": 0,
    "overwrite": False,
    "output_compression": None,
}
meta = {"name": "harmony_integrate"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from anndata_io import read_modality, write_modality

logger = setup_logger()

dat = read_modality(par["input"], par["modality"], logger)

logger.info(par)

if par["obsm_input"] not in dat.obsm.keys():
    raise ValueError(
        f"obsm slot '{par['obsm_input']}' not found in modality '{par['modality']}'."
    )

for covariate in par["obs_covariates"]:
    if covariate not in dat.obs.columns:
        raise ValueError(
            f"obs column '{covariate}' not found in modality '{par['modality']}'."
        )

if par["obsm_output"] in dat.obsm and not par["overwrite"]:
    raise ValueError(
        f"Requested to create field {par['obsm_output']} in .obsm for modality "
        f"{par['modality']}, but field already exists."
    )

# rsc.pp.harmony_integrate reads the embedding from .obsm[basis] and writes the
# corrected embedding to .obsm[adjusted_basis]. It does not transform the main
# matrix, so the obsm representation is transferred to GPU internally - call it
# directly, no on_gpu pre-stage required.

logger.info("Integrating with Harmony.")
rsc.pp.harmony_integrate(
    dat,
    key=par["obs_covariates"],
    basis=par["obsm_input"],
    adjusted_basis=par["obsm_output"],
    theta=par["theta"],
    n_clusters=par["n_clusters"],
    max_iter_harmony=par["max_iter_harmony"],
    random_state=par["random_state"],
)

write_modality(
    dat, par["output"], par["input"], par["modality"], par["output_compression"], logger
)
