import sys
import rapids_singlecell as rsc

## VIASH START
par = {
    "input": "resources_test/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu",
    "output": "output.h5mu",
    "modality": "rna",
    "input_layer": None,
    "var_name_filter": "highly_variable",
    "varm_name": "hvg",
    "flavor": "seurat",
    "n_top_features": None,
    "min_mean": 0.0125,
    "max_mean": 3.0,
    "min_disp": 0.5,
    "max_disp": None,
    "span": 0.3,
    "n_bins": 20,
    "theta": 100,
    "clip": None,
    "chunksize": 1000,
    "n_samples": 10000,
    "obs_batch_key": None,
    "check_values": True,
}
meta = {"name": "highly_variable_genes"}
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

# Flavors that require n_top_features
flavors_requiring_n_top = {
    "seurat_v3",
    "seurat_v3_paper",
    "pearson_residuals",
    "poisson_gene_selection",
}
if par["flavor"] in flavors_requiring_n_top and not par["n_top_features"]:
    raise ValueError(
        f"When flavor is set to '{par['flavor']}', you are required to set "
        "'n_top_features'."
    )

# rapids_singlecell.pp.highly_variable_genes writes to adata.var / adata.uns
# rather than transforming a matrix, so we do not pre-stage an output layer;
# we just point it at the requested input matrix via `layer`.

# Build kwargs, only forwarding optional values that are set so we keep the
# upstream defaults for everything else.
hvg_kwargs = {
    "layer": par["input_layer"],
    "flavor": par["flavor"],
    "min_mean": par["min_mean"],
    "max_mean": par["max_mean"],
    "min_disp": par["min_disp"],
    "span": par["span"],
    "n_bins": par["n_bins"],
    "theta": par["theta"],
    "chunksize": par["chunksize"],
    "n_samples": par["n_samples"],
    "check_values": par["check_values"],
}
if par["n_top_features"] is not None:
    hvg_kwargs["n_top_genes"] = par["n_top_features"]
if par["max_disp"] is not None:
    hvg_kwargs["max_disp"] = par["max_disp"]
if par["clip"] is not None:
    hvg_kwargs["clip"] = par["clip"]
if par["obs_batch_key"] is not None:
    hvg_kwargs["batch_key"] = par["obs_batch_key"]

with on_gpu(dat, logger, layer=par["input_layer"]):
    logger.info("Computing highly variable genes.")
    rsc.pp.highly_variable_genes(dat, **hvg_kwargs)

# rapids-singlecell writes a fixed set of columns to .var. Rename/move them so
# the output matches the configured --var_name_filter and --varm_name layout.
hvg_var_columns = [
    "highly_variable",
    "means",
    "dispersions",
    "dispersions_norm",
    "variances",
    "variances_norm",
    "residual_variances",
    "highly_variable_rank",
    "highly_variable_nbatches",
    "highly_variable_intersection",
]
present_columns = [c for c in hvg_var_columns if c in dat.var.columns]

if par["varm_name"] and present_columns:
    # Copy the full per-gene HVG output into .varm, including the boolean
    # "highly_variable" flag.
    dat.varm[par["varm_name"]] = dat.var[present_columns].copy()

if par["var_name_filter"] and par["var_name_filter"] != "highly_variable":
    dat.var[par["var_name_filter"]] = dat.var["highly_variable"]

if par["varm_name"] and present_columns:
    # .var only keeps the boolean filter column; the metrics live in .varm.
    dat.var = dat.var.drop(
        columns=[c for c in present_columns if c != par["var_name_filter"]]
    )

write_modality(
    dat, par["output"], par["input"], par["modality"], par["output_compression"], logger
)
