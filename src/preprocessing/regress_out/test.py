import sys
import subprocess
import pytest
import mudata as mu
import numpy as np

## VIASH START
meta = {
    "name": "regress_out",
    "resources_dir": "resources_test/",
    "config": "src/preprocessing/regress_out/config.vsh.yaml",
    "executable": "target/docker/preprocessing/regress_out/regress_out",
}
## VIASH END

input = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"


@pytest.fixture
def clean_input(random_h5mu_path):
    """Ensure the test data carries a numeric .obs covariate to regress out.

    The "_mms" input normally has QC columns, but to be safe we add a
    numeric total_counts column when it is missing, then write the data to a
    fresh path for the component to consume."""
    mu_in = mu.read_h5mu(input)
    rna = mu_in.mod["rna"]
    if "total_counts" not in rna.obs.columns:
        rna.obs["total_counts"] = np.asarray(rna.X.sum(axis=1)).ravel()
    path = random_h5mu_path()
    mu_in.write(path)
    return path


def test_run(run_component, random_h5mu_path, clean_input):
    output = random_h5mu_path()
    run_component(
        [
            "--input",
            clean_input,
            "--output",
            output,
            "--output_compression",
            "gzip",
            "--obs_keys",
            "total_counts",
        ]
    )

    assert output.is_file(), "No output was created."

    mu_input = mu.read_h5mu(clean_input)
    mu_output = mu.read_h5mu(output)

    assert "rna" in mu_output.mod, 'Output should contain data.mod["rna"].'

    rna_in = mu_input.mod["rna"]
    rna_out = mu_output.mod["rna"]

    assert rna_in.shape == rna_out.shape, "Should have same shape as before"
    assert not np.allclose(
        np.asarray(rna_in.X.todense() if hasattr(rna_in.X, "todense") else rna_in.X),
        np.asarray(rna_out.X.todense() if hasattr(rna_out.X, "todense") else rna_out.X),
    ), "Expression should have changed"


def test_output_layer(run_component, random_h5mu_path, clean_input):
    """Writing the result to a named output layer should leave X untouched
    and store the regressed matrix in that layer."""
    output = random_h5mu_path()
    run_component(
        [
            "--input",
            clean_input,
            "--output",
            output,
            "--output_compression",
            "gzip",
            "--obs_keys",
            "total_counts",
            "--output_layer",
            "regressed",
        ]
    )

    rna_in = mu.read_h5mu(clean_input).mod["rna"]
    rna_out = mu.read_h5mu(output).mod["rna"]

    in_x = np.asarray(rna_in.X.todense() if hasattr(rna_in.X, "todense") else rna_in.X)
    out_x = np.asarray(
        rna_out.X.todense() if hasattr(rna_out.X, "todense") else rna_out.X
    )
    assert np.allclose(out_x, in_x), (
        "X should be untouched when writing to an output layer"
    )
    assert "regressed" in rna_out.layers, "Output layer should be present"
    assert rna_out.layers["regressed"].shape == rna_out.shape


def test_missing_obs_key_raises(run_component, random_h5mu_path, clean_input):
    output = random_h5mu_path()
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                clean_input,
                "--output",
                output,
                "--output_compression",
                "gzip",
                "--obs_keys",
                "does_not_exist",
            ]
        )
    assert "does_not_exist" in err.value.stdout.decode("utf-8")


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
