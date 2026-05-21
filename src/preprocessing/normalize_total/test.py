import sys
import pytest
import mudata as mu
import numpy as np
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


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
