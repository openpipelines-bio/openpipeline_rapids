import sys
import subprocess
import pytest
import mudata as mu
import numpy as np

## VIASH START
meta = {
    "name": "scale",
    "resources_dir": "resources_test/",
    "config": "src/preprocessing/scale/config.vsh.yaml",
    "executable": "target/docker/preprocessing/scale/scale",
}
## VIASH END

input = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"


def _to_dense(matrix):
    return matrix.toarray() if hasattr(matrix, "toarray") else np.asarray(matrix)


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
    assert np.mean(_to_dense(rna_in.X)) != np.mean(_to_dense(rna_out.X)), (
        "Expression should have changed"
    )


def test_output_layer(run_component, random_h5mu_path):
    """Writing the result to a named output layer should leave X untouched
    and store the scaled matrix in that layer."""
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
            "scaled",
        ]
    )

    rna_in = mu.read_h5mu(input).mod["rna"]
    rna_out = mu.read_h5mu(output).mod["rna"]

    assert "scaled" in rna_out.layers, "Output layer should have been created."
    assert np.allclose(_to_dense(rna_out.X), _to_dense(rna_in.X)), (
        "X should be untouched when writing to an output layer"
    )
    assert not np.allclose(_to_dense(rna_out.layers["scaled"]), _to_dense(rna_in.X)), (
        "Output layer should hold the scaled matrix"
    )


def test_max_value(run_component, random_h5mu_path):
    """max_value should clip the scaled values so none exceed the cap."""
    max_value = 1.0
    output = random_h5mu_path()
    run_component(
        [
            "--input",
            input,
            "--output",
            output,
            "--output_compression",
            "gzip",
            "--max_value",
            str(max_value),
        ]
    )

    rna_out = mu.read_h5mu(output).mod["rna"]
    assert np.max(_to_dense(rna_out.X)) <= max_value + 1e-5, (
        "Scaled values should be clipped at max_value"
    )


def test_raise_if_input_layer_missing(run_component, random_h5mu_path):
    """A missing --input_layer should make the component fail."""
    output = random_h5mu_path()
    with pytest.raises(subprocess.CalledProcessError):
        run_component(
            [
                "--input",
                input,
                "--output",
                output,
                "--output_compression",
                "gzip",
                "--input_layer",
                "does_not_exist",
            ]
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


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
