import mudata as mu
from compress_h5mu import write_h5ad_to_h5mu_with_compression


def read_modality(par, logger):
    logger.info("Reading modality %s from %s", par["modality"], par["input"])
    dat = mu.read_h5ad(par["input"], mod=par["modality"])
    assert dat.var_names.is_unique, (
        "The var_names of the input modality must be unique."
    )
    return dat


def write_modality(par, dat, logger):
    logger.info(
        "Writing to file to %s with compression %s",
        par["output"],
        par["output_compression"],
    )
    write_h5ad_to_h5mu_with_compression(
        par["output"], par["input"], par["modality"], dat, par["output_compression"]
    )
