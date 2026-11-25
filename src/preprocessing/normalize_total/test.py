import logging
import sys
import pytest
import mudata as mu
import numpy as np
import scipy.sparse
from openpipeline_testutils.asserters import assert_annotation_objects_equal

## VIASH START
meta = {
    "name": "normalize_total",
    "resources_dir": "resources_test/",
    "config": "src/preprocessing/normalize_total/config.vsh.yaml",
    "executable": "target/docker/preprocessing/normalize_total/normalize_total",
}
## VIASH END

input = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"


def test_run(run_component, random_h5mu_path):
    output = random_h5mu_path()
    cmd_pars = [
        "--input",
        input,
        "--output",
        output,
        "--output_compression",
        "gzip",
    ]
    run_component(cmd_pars)

    assert output.is_file(), "No output was created."

    mu_input = mu.read_h5mu(input)
    mu_output = mu.read_h5mu(output)

    assert "rna" in mu_output.mod, 'Output should contain data.mod["rna"].'

    rna_in = mu_input.mod["rna"]
    rna_out = mu_output.mod["rna"]

    assert rna_in.shape == rna_out.shape, "Should have same shape as before"
    assert np.mean(rna_in.X) != np.mean(rna_out.X), "Expression should have changed"

    nz_row, nz_col = rna_in.X.nonzero()
    row_corr = np.corrcoef(
        rna_in.X[nz_row[0], :].toarray().flatten(),
        rna_out.X[nz_row[0], :].toarray().flatten(),
    )[0, 1]
    col_corr = np.corrcoef(
        rna_in.X[:, nz_col[0]].toarray().flatten(),
        rna_out.X[:, nz_col[0]].toarray().flatten(),
    )[0, 1]
    assert row_corr > 0.1
    assert col_corr > 0.1

    # Copy over X so that the rest of the objects can be compared
    mu_input["rna"].X = mu_output["rna"].X
    assert_annotation_objects_equal(mu_input, mu_output)


def test_target_sum(run_component, random_h5mu_path):
    output = random_h5mu_path()
    cmd_pars = [
        "--input",
        input,
        "--output",
        output,
        "--output_compression",
        "gzip",
        "--target_sum",
        "10000",
    ]
    run_component(cmd_pars)

    assert output.is_file(), "No output was created."

    mu_output = mu.read_h5mu(output)
    mu_input = mu.read_h5mu(input)
    rna_out = mu_output.mod["rna"]

    assert np.all(np.abs(rna_out.X.sum(axis=1) - 10000) < 1), (
        "Expression should have changed"
    )

    # Copy over X so that the rest of the object can be compared
    mu_input["rna"].X = mu_output["rna"].X
    assert_annotation_objects_equal(mu_input, mu_output)


def test_exclude_highly_expressed(run_component, random_h5mu_path):
    """With an aggressive max_fraction, excluding highly expressed genes
    should produce a different result than the default normalization."""
    output_default = random_h5mu_path()
    output_excl = random_h5mu_path()

    run_component(
        [
            "--input",
            input,
            "--output",
            output_default,
            "--output_compression",
            "gzip",
            "--target_sum",
            "10000",
        ]
    )
    run_component(
        [
            "--input",
            input,
            "--output",
            output_excl,
            "--output_compression",
            "gzip",
            "--target_sum",
            "10000",
            "--exclude_highly_expressed",
            "--max_fraction",
            "0.001",
        ]
    )

    rna_default = mu.read_h5mu(output_default).mod["rna"]
    rna_excl = mu.read_h5mu(output_excl).mod["rna"]

    assert rna_default.shape == rna_excl.shape
    assert not np.allclose(rna_default.X.toarray(), rna_excl.X.toarray()), (
        "exclude_highly_expressed should change the normalization result"
    )


def test_input_layer(run_component, random_h5mu_path):
    """Normalizing a named layer should leave X untouched and write the
    normalized matrix back into that layer."""
    mu_orig = mu.read_h5mu(input)
    mu_orig.mod["rna"].layers["counts"] = mu_orig.mod["rna"].X.copy()
    input_with_layer = random_h5mu_path()
    mu_orig.write(input_with_layer)

    output = random_h5mu_path()
    run_component(
        [
            "--input",
            str(input_with_layer),
            "--output",
            output,
            "--output_compression",
            "gzip",
            "--target_sum",
            "10000",
            "--input_layer",
            "counts",
        ]
    )

    rna_out = mu.read_h5mu(output).mod["rna"]
    rna_orig = mu_orig.mod["rna"]

    assert np.allclose(rna_out.X.toarray(), rna_orig.X.toarray()), (
        "X should be untouched when normalizing a layer"
    )
    layer_sums = rna_out.layers["counts"].sum(axis=1)
    assert np.all(np.abs(layer_sums - 10000) < 1), (
        "Layer sum per cell should equal target_sum"
    )


def test_output_layer(run_component, random_h5mu_path):
    """Writing the result to a named output layer should leave X untouched
    and store the normalized matrix in that layer."""
    output = random_h5mu_path()
    run_component(
        [
            "--input",
            input,
            "--output",
            output,
            "--output_compression",
            "gzip",
            "--target_sum",
            "10000",
            "--output_layer",
            "normalized",
        ]
    )

    rna_in = mu.read_h5mu(input).mod["rna"]
    rna_out = mu.read_h5mu(output).mod["rna"]

    assert np.allclose(rna_out.X.toarray(), rna_in.X.toarray()), (
        "X should be untouched when writing to an output layer"
    )
    layer_sums = rna_out.layers["normalized"].sum(axis=1)
    assert np.all(np.abs(layer_sums - 10000) < 1), (
        "Output layer sum per cell should equal target_sum"
    )


def test_empty_input_layer_uses_x(run_component, random_h5mu_path):
    """An empty --input_layer means ".X"; it must not be treated as a layer name."""
    output = random_h5mu_path()
    run_component(
        [
            "--input",
            input,
            "--output",
            output,
            "--input_layer",
            "",
        ]
    )

    assert output.is_file(), "No output was created."
    rna_in = mu.read_h5mu(input).mod["rna"]
    rna_out = mu.read_h5mu(output).mod["rna"]
    assert rna_in.shape == rna_out.shape
    assert np.mean(rna_in.X) != np.mean(rna_out.X), "Expression should have changed"


def test_on_gpu_restores_cpu_matrices_on_error():
    """When the GPU block fails, on_gpu puts the untouched CPU matrices back so
    that a retry starts from the original input."""
    sys.path.append(meta["resources_dir"])
    import rapids_singlecell as rsc
    from gpu import on_gpu

    rna = mu.read_h5mu(input).mod["rna"]
    rna.layers["counts"] = rna.X.copy()
    original_x = rna.X.copy()

    with pytest.raises(RuntimeError, match="Step failed on the GPU"):
        with on_gpu(rna, convert_all=True):
            rsc.pp.normalize_total(rna, target_sum=10000)
            raise RuntimeError("Step failed on the GPU")

    for name, matrix in [("X", rna.X), ("counts", rna.layers["counts"])]:
        assert scipy.sparse.issparse(matrix), f"{name} should be back on the CPU"
        assert (matrix != original_x).nnz == 0, f"{name} should be unchanged"


def test_host_memory_error_is_not_gpu_out_of_memory():
    """A host allocation failure must not trigger the managed-memory retry."""
    sys.path.append(meta["resources_dir"])
    from gpu import _is_gpu_out_of_memory

    with pytest.raises(MemoryError) as excinfo:
        np.empty(2**60, dtype=np.uint8)
    assert not _is_gpu_out_of_memory(excinfo.value)


def test_managed_memory_fallback(caplog):
    """A step that needs more memory than the GPU has free runs out of memory,
    and succeeds when it is retried with RMM managed memory."""
    sys.path.append(meta["resources_dir"])
    import cupy as cp
    import rmm
    from gpu import run_with_managed_memory_fallback

    free_bytes, _ = cp.cuda.runtime.memGetInfo()
    oversized_bytes = int(free_bytes * 1.5)
    memory_resources = []

    def allocate_more_than_free_vram():
        memory_resources.append(rmm.mr.get_current_device_resource())
        # cp.empty does not touch the pages, so the oversubscribed managed
        # allocation does not need host RAM to back it.
        return cp.empty(oversized_bytes, dtype=cp.uint8).nbytes

    logger = logging.getLogger("test_managed_memory_fallback")
    with caplog.at_level(logging.WARNING):
        allocated = run_with_managed_memory_fallback(
            allocate_more_than_free_vram, logger
        )

    assert allocated == oversized_bytes
    assert len(memory_resources) == 2, "The step should run twice"
    assert not isinstance(memory_resources[0], rmm.mr.ManagedMemoryResource)
    assert isinstance(memory_resources[1], rmm.mr.ManagedMemoryResource)
    assert "Retrying with RMM managed memory" in caplog.text


def test_on_gpu_releases_cpu_matrices_under_managed_memory():
    """Once managed memory is enabled, on_gpu keeps no CPU copy of the input, so
    host RAM is left for the data managed memory pages out of VRAM.

    Enables managed memory for the rest of the test process, so it runs last.
    """
    sys.path.append(meta["resources_dir"])
    import gc
    import weakref
    import gpu

    rna = mu.read_h5mu(input).mod["rna"]
    rna.X = rna.X.astype(np.float32)
    cpu_x = weakref.ref(rna.X)

    if not gpu._managed_memory_enabled:
        gpu._enable_managed_memory()
    with gpu.on_gpu(rna):
        gc.collect()
        assert cpu_x() is None, "The CPU copy of X should have been released"

    assert scipy.sparse.issparse(rna.X), "X should be back on the CPU"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
