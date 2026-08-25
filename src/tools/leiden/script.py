import sys
import rapids_singlecell as rsc
import pandas as pd

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "obsp_connectivities": "connectivities",
    "obsm_name": "leiden",
    "resolution": [1.0, 0.25],
    "n_iterations": 100,
    "theta": 1.0,
    "random_state": 0,
    "use_weights": True,
    "output_compression": None,
}
meta = {"name": "leiden"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from anndata_io import read_modality, write_modality

logger = setup_logger()

dat = read_modality(par["input"], par["modality"], logger)

logger.info(par)

if par["obsp_connectivities"] not in dat.obsp.keys():
    raise ValueError(
        f"Could not find .obsp key '{par['obsp_connectivities']}' in modality "
        f"'{par['modality']}'. Run a neighbors component first."
    )

# rapids-singlecell's tl.leiden writes cluster labels to .obs and parameters
# to .uns. It does not transform a matrix, so the GPU transfer for X/layer is
# not needed; cuGraph builds its own graph from the adjacency matrix.

# Run Leiden once per resolution and collect the labels into a DataFrame so
# that the output matches the openpipeline CPU leiden component interface
# (.obsm[obsm_name] with one column per resolution).
results = {}
for resolution in par["resolution"]:
    logger.info("Running Leiden clustering at resolution %s.", resolution)
    rsc.tl.leiden(
        dat,
        resolution=resolution,
        random_state=par["random_state"],
        theta=par["theta"],
        n_iterations=par["n_iterations"],
        use_weights=par["use_weights"],
        obsp=par["obsp_connectivities"],
        key_added=par["obsm_name"],
    )
    results[str(resolution)] = dat.obs[par["obsm_name"]].copy()
    # Drop the per-resolution scratch column so it doesn't leak into output.obs
    del dat.obs[par["obsm_name"]]

# Also drop the .uns entry that rsc.tl.leiden writes - we keep results in
# .obsm to match the openpipeline interface.
if par["obsm_name"] in dat.uns:
    del dat.uns[par["obsm_name"]]

dat.obsm[par["obsm_name"]] = pd.DataFrame(results, index=dat.obs_names)

write_modality(
    dat, par["output"], par["input"], par["modality"], par["output_compression"], logger
)
