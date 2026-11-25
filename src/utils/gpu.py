import gc
from contextlib import contextmanager
import cupy as cp
import numpy as np
import rapids_singlecell as rsc
import rmm
from rmm.allocators.cupy import rmm_cupy_allocator

# Substrings (lowercase) that identify out-of-memory errors raised by the CUDA
# libraries rapids-singlecell builds on: RMM and cuML/cuGraph/cuVS raise
# ``std::bad_alloc: out_of_memory: CUDA error ...``, cupy raises
# ``cudaErrorMemoryAllocation: out of memory`` and cuSPARSE/cuBLAS report
# ``*_STATUS_ALLOC_FAILED``. A bare ``bad_alloc`` is not on the list: host-side
# C++ allocation failures report it too.
_GPU_OOM_MARKERS = (
    "out of memory",
    "out_of_memory",
    "alloc_failed",
    "memoryallocation",
)

_managed_memory_enabled = False


def _ensure_gpu_compatible_dtype(matrix):
    """Cast integer matrices to float32 so they can be moved to the GPU.

    cupy sparse matrices only support bool/float/complex dtypes, so integer
    count matrices (e.g. raw counts) cannot be transferred as-is. The
    rapids-singlecell operations that run on the transferred data produce float
    values anyway, so casting integer counts to float loses no information.
    """
    if matrix is not None and np.issubdtype(matrix.dtype, np.integer):
        return matrix.astype("float32")
    return matrix


def _cast_transfer_targets(adata, kwargs):
    """Cast the matrices ``anndata_to_GPU`` will transfer to a GPU-safe dtype.

    Mirrors what ``anndata_to_GPU`` moves for the given keyword arguments so
    that matrices it does not touch keep their original dtype.
    """
    if kwargs.get("convert_all"):
        adata.X = _ensure_gpu_compatible_dtype(adata.X)
        for key in list(adata.layers.keys()):
            adata.layers[key] = _ensure_gpu_compatible_dtype(adata.layers[key])
    elif kwargs.get("layer"):
        layer = kwargs["layer"]
        adata.layers[layer] = _ensure_gpu_compatible_dtype(adata.layers[layer])
    else:
        adata.X = _ensure_gpu_compatible_dtype(adata.X)


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

    Integer matrices among the transferred ``.X``/``.layers`` are cast to
    float32 first, since cupy sparse cannot represent integer dtypes.

    If the transfer to the GPU, the block or the transfer back raises, the CPU
    matrices that were moved are put back before the exception propagates, so
    ``adata`` is left in its pre-transfer state and holds no GPU memory. This
    lets ``run_with_managed_memory_fallback`` rerun the block from clean input.
    """
    # Components treat an empty --layer as "use .X", and pass it through
    # explicitly rather than omitting the argument. An empty string is not a
    # layer name, so drop it instead of forwarding it to anndata_to_GPU.
    if "layer" in kwargs and not kwargs["layer"]:
        del kwargs["layer"]

    slot_keys = {
        attr: [keys] if isinstance(keys, str) else list(keys)
        for attr, keys in (slots or {}).items()
    }
    transfer_adata = bool(kwargs) or not slot_keys

    def move_slots(transfer_func):
        for attr, keys in slot_keys.items():
            mapping = getattr(adata, attr)
            for key in keys:
                mapping[key] = transfer_func(mapping[key])

    description = _describe_transfer(transfer_adata, kwargs, slot_keys)

    if logger is not None:
        logger.info("Transferring %s to GPU.", description)
    if transfer_adata:
        _cast_transfer_targets(adata, kwargs)
    cpu_originals = _collect_transfer_targets(adata, transfer_adata, kwargs, slot_keys)

    try:
        if transfer_adata:
            rsc.get.anndata_to_GPU(adata, **kwargs)
        move_slots(rsc.get.X_to_GPU)

        yield adata

        if logger is not None:
            logger.info("Transferring %s back to CPU.", description)
        if transfer_adata:
            rsc.get.anndata_to_CPU(adata)
        move_slots(rsc.get.X_to_CPU)
    except BaseException:
        _restore_transfer_targets(adata, cpu_originals)
        raise


def _collect_transfer_targets(adata, transfer_adata, kwargs, slot_keys):
    """References to the CPU objects ``on_gpu`` is about to move to the GPU.

    Returns a list of ``(attr, key, value)`` tuples, where ``key`` is ``None``
    for ``.X``. Mirrors what ``anndata_to_GPU`` and ``X_to_GPU`` replace for
    the given arguments.
    """
    targets = []
    if transfer_adata:
        if kwargs.get("convert_all"):
            targets.append(("X", None, adata.X))
            # anndata >= 0.13 also lists .X under the None key of .layers; it is
            # already recorded above.
            targets.extend(
                ("layers", key, adata.layers[key])
                for key in adata.layers
                if key is not None
            )
        elif kwargs.get("layer"):
            layer = kwargs["layer"]
            targets.append(("layers", layer, adata.layers[layer]))
        else:
            targets.append(("X", None, adata.X))
    for attr, keys in slot_keys.items():
        mapping = getattr(adata, attr)
        targets.extend((attr, key, mapping[key]) for key in keys)
    return targets


def _restore_transfer_targets(adata, targets):
    """Put back the CPU objects collected by ``_collect_transfer_targets``."""
    for attr, key, value in targets:
        if key is None:
            setattr(adata, attr, value)
        else:
            getattr(adata, attr)[key] = value


def _is_gpu_out_of_memory(exc):
    """Whether ``exc`` signals that the GPU ran out of memory.

    cupy raises ``OutOfMemoryError``; RMM and the CUDA libraries raise
    ``MemoryError`` or ``RuntimeError`` with a CUDA allocation message. A
    ``MemoryError`` without such a message, such as numpy's "Unable to
    allocate ...", is a host allocation failure that managed memory cannot fix.
    """
    if isinstance(exc, cp.cuda.memory.OutOfMemoryError):
        return True
    message = str(exc).lower()
    return any(marker in message for marker in _GPU_OOM_MARKERS)


def _enable_managed_memory():
    """Switch RMM (and cupy, through RMM) to CUDA managed memory.

    Managed memory oversubscribes the GPU by paging data to host RAM, so
    datasets larger than VRAM can be processed at the cost of speed. RMM's
    current device resource is also what cuML, cuGraph and cuVS allocate from,
    so their allocations are covered as well. See
    https://rapids-singlecell.scverse.org/en/latest/memory_management.html
    """
    global _managed_memory_enabled
    rmm.reinitialize(
        managed_memory=True,
        pool_allocator=False,
        devices=cp.cuda.Device().id,
    )
    cp.cuda.set_allocator(rmm_cupy_allocator)
    _managed_memory_enabled = True


def run_with_managed_memory_fallback(step, logger):
    """Run ``step()``; if the GPU runs out of memory, retry once with managed memory.

    ``step`` must leave its input data unchanged when it fails, so that it can
    be rerun from the start. Steps that move data with ``on_gpu`` get this from
    ``on_gpu``, which restores the CPU matrices on failure. Steps that call a
    rapids-singlecell function on CPU data get it from rapids-singlecell, which
    transfers copies to the GPU internally.

    Once managed memory is enabled it stays enabled for the rest of the
    process, and later out-of-memory errors propagate.
    """
    try:
        return step()
    except Exception as exc:
        if _managed_memory_enabled or not _is_gpu_out_of_memory(exc):
            raise
        logger.warning(
            "GPU ran out of memory (%s). Retrying with RMM managed memory, which "
            "spills to host RAM and is slower.",
            exc,
        )
    # Reinitializing RMM while earlier allocations are alive is undefined
    # behaviour. Retrying outside the except block releases the traceback and
    # the GPU buffers its frames still reference; collect them before switching.
    gc.collect()
    cp.get_default_memory_pool().free_all_blocks()
    _enable_managed_memory()
    return step()


def _describe_transfer(transfer_adata, kwargs, slot_keys):
    """Human-readable summary of what ``on_gpu`` moves, for logging."""
    parts = []
    if transfer_adata:
        if kwargs.get("convert_all"):
            parts.append(".X and all layers")
        elif "layer" in kwargs:
            parts.append(f"layer {kwargs['layer']!r}")
        else:
            parts.append(".X and .layers")
    for attr, keys in slot_keys.items():
        keys_str = ", ".join(repr(key) for key in keys)
        parts.append(f".{attr}[{keys_str}]")
    return ", ".join(parts) if parts else "nothing"
