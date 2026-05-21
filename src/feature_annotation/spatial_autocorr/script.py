import sys
import rapids_singlecell as rsc
import mudata as mu

## VIASH START
par = {
    "input": "resources_test/xenium/xenium_tiny.qc.neighbors.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "layer": None,
    "obsp_connectivities": "spatial_connectivities",
    "genes": None,
    "use_raw": False,
    "mode": "moran",
    "n_perms": None,
    "transformation": True,
    "two_tailed": False,
    "corr_method": "fdr_bh",
    "use_sparse": True,
    "output_compression": None,
}
meta = {"name": "spatial_autocorr"}
## VIASH END

sys.path.append(meta["resources_dir"])
from setup_logger import setup_logger
from compress_h5mu import write_h5ad_to_h5mu_with_compression

logger = setup_logger()

logger.info("Reading modality %s from %s", par["modality"], par["input"])
dat = mu.read_h5ad(par["input"], mod=par["modality"])

logger.info(par)

if par["obsp_connectivities"] not in dat.obsp.keys():
    raise ValueError(
        f"Could not find .obsp key '{par['obsp_connectivities']}' in modality "
        f"'{par['modality']}'. Run a spatial neighbors component first."
    )

if par["layer"] and par["layer"] not in dat.layers.keys():
    raise ValueError(
        f"Layer '{par['layer']}' not found in modality '{par['modality']}'."
    )

genes = par["genes"]
if genes is not None and len(genes) == 0:
    genes = None

# An empty --corr_method disables multiple-testing correction.
corr_method = par["corr_method"] or None

logger.info("Computing spatial autocorrelation (%s).", par["mode"])
rsc.gr.spatial_autocorr(
    dat,
    connectivity_key=par["obsp_connectivities"],
    genes=genes,
    mode=par["mode"],
    transformation=par["transformation"],
    n_perms=par["n_perms"],
    two_tailed=par["two_tailed"],
    corr_method=corr_method,
    layer=par["layer"],
    use_raw=par["use_raw"],
    use_sparse=par["use_sparse"],
)

result_key = "moranI" if par["mode"] == "moran" else "gearyC"
if result_key in dat.uns:
    df = dat.uns[result_key]
    if not df.empty:
        stat_col = "I" if par["mode"] == "moran" else "C"
        logger.info("Top spatially variable features:")
        logger.info(df.sort_values(by=stat_col, ascending=False).head().to_string())
else:
    logger.warning("Expected key '%s' not found in .uns after computation.", result_key)

logger.info(
    "Writing to file to %s with compression %s",
    par["output"],
    par["output_compression"],
)
write_h5ad_to_h5mu_with_compression(
    par["output"], par["input"], par["modality"], dat, par["output_compression"]
)
