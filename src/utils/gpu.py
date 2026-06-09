from contextlib import contextmanager
import rapids_singlecell as rsc


@contextmanager
def on_gpu(adata, logger=None, *, slots=None, **kwargs):
    """Move an AnnData to the GPU for the duration of the block, then back.

    Extra keyword arguments are forwarded to ``rsc.get.anndata_to_GPU``
    (e.g. ``layer=...`` or ``convert_all=True``).

    ``anndata_to_GPU`` only transfers ``.X`` and ``.layers``; it never touches
    aligned mappings such as ``.obsm`` or ``.varm``. Pass ``slots`` to move
    individual entries of those via ``X_to_GPU``, restoring them with
    ``X_to_CPU`` on exit. ``slots`` maps an AnnData attribute name to the
    key(s) within it, e.g. ``{"obsm": "X_pca"}`` or
    ``{"obsm": ["X_pca", "X_umap"], "varm": ["loadings"]}``. When only ``slots``
    is requested, ``.X``/``.layers`` are left on the CPU.
    """
    slot_keys = {
        attr: [keys] if isinstance(keys, str) else list(keys)
        for attr, keys in (slots or {}).items()
    }
    transfer_adata = bool(kwargs) or not slot_keys

    def move_slots(transfer):
        for attr, keys in slot_keys.items():
            mapping = getattr(adata, attr)
            for key in keys:
                mapping[key] = transfer(mapping[key])

    if logger is not None:
        logger.info("Transferring data to GPU.")
    if transfer_adata:
        rsc.get.anndata_to_GPU(adata, **kwargs)
    move_slots(rsc.get.X_to_GPU)

    yield adata

    if logger is not None:
        logger.info("Transferring data back to CPU.")
    if transfer_adata:
        rsc.get.anndata_to_CPU(adata)
    move_slots(rsc.get.X_to_CPU)
