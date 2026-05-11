import sys
import pytest
import mudata as mu
import numpy as np
from openpipeline_testutils.asserters import assert_annotation_objects_equal

## VIASH START
meta = {
    "name": "log1p",
    "resources_dir": "resources_test/",
    "config": "src/preprocessing/log1p/config.vsh.yaml",
    "executable": "target/docker/preprocessing/log1p/log1p",
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

    # log1p(0) == 0, so zero entries must remain zero
    in_dense = rna_in.X.toarray()
    out_dense = rna_out.X.toarray()
    assert np.all(out_dense[in_dense == 0] == 0), (
        "Zero entries should remain zero after log1p"
    )
    # Non-zero entries should equal natural log(1 + x)
    nz_mask = in_dense != 0
    assert np.allclose(out_dense[nz_mask], np.log1p(in_dense[nz_mask])), (
        "Non-zero entries should equal log1p of the input"
    )

    # Copy over X so that the rest of the objects can be compared
    mu_input["rna"].X = mu_output["rna"].X
    assert_annotation_objects_equal(mu_input, mu_output)


def test_base(run_component, random_h5mu_path):
    """Using a custom base should rescale the result by 1 / log(base)."""
    output_natural = random_h5mu_path()
    output_base10 = random_h5mu_path()

    run_component(
        [
            "--input",
            input,
            "--output",
            output_natural,
            "--output_compression",
            "gzip",
        ]
    )
    run_component(
        [
            "--input",
            input,
            "--output",
            output_base10,
            "--output_compression",
            "gzip",
            "--base",
            "10",
        ]
    )

    natural = mu.read_h5mu(output_natural).mod["rna"].X.toarray()
    base10 = mu.read_h5mu(output_base10).mod["rna"].X.toarray()

    nz_mask = natural != 0
    assert np.allclose(base10[nz_mask], natural[nz_mask] / np.log(10)), (
        "base=10 result should equal natural log result divided by log(10)"
    )


def test_input_layer(run_component, random_h5mu_path):
    """Log-transforming a named layer should leave X untouched and write the
    transformed matrix back into that layer."""
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
            "--input_layer",
            "counts",
        ]
    )

    rna_out = mu.read_h5mu(output).mod["rna"]
    rna_orig = mu_orig.mod["rna"]

    assert np.allclose(rna_out.X.toarray(), rna_orig.X.toarray()), (
        "X should be untouched when transforming a layer"
    )
    in_dense = rna_orig.layers["counts"].toarray()
    out_dense = rna_out.layers["counts"].toarray()
    assert np.allclose(out_dense, np.log1p(in_dense)), (
        "Layer should equal log1p of the input counts"
    )


def test_output_layer(run_component, random_h5mu_path):
    """Writing the result to a named output layer should leave X untouched
    and store the log-transformed matrix in that layer."""
    output = random_h5mu_path()
    run_component(
        [
            "--input",
            input,
            "--output",
            output,
            "--output_compression",
            "gzip",
            "--output_layer",
            "log_normalized",
        ]
    )

    rna_in = mu.read_h5mu(input).mod["rna"]
    rna_out = mu.read_h5mu(output).mod["rna"]

    assert np.allclose(rna_out.X.toarray(), rna_in.X.toarray()), (
        "X should be untouched when writing to an output layer"
    )
    in_dense = rna_in.X.toarray()
    out_dense = rna_out.layers["log_normalized"].toarray()
    assert np.allclose(out_dense, np.log1p(in_dense)), (
        "Output layer should equal log1p of the input X"
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
