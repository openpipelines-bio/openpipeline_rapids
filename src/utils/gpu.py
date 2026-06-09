from contextlib import contextmanager
import rapids_singlecell as rsc


@contextmanager
def on_gpu(adata, logger=None, **kwargs):
    """Move an AnnData to the GPU for the duration of the block, then back.

    Extra keyword arguments are forwarded to ``rsc.get.anndata_to_GPU``
    (e.g. ``layer=...`` or ``convert_all=True``).
    """
    if logger is not None:
        logger.info("Transferring data to GPU.")
    rsc.get.anndata_to_GPU(adata, **kwargs)
    yield adata
    if logger is not None:
        logger.info("Transferring data back to CPU.")
    rsc.get.anndata_to_CPU(adata)
