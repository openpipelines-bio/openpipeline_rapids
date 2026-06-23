import sys
import pytest
import subprocess
import mudata as mu

## VIASH START
meta = {
    "name": "tsne",
    "resources_dir": "resources_test/",
    "config": "src/tools/tsne/config.vsh.yaml",
    "executable": "target/docker/tools/tsne/tsne",
}
## VIASH END

input = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"


def test_run(run_component, random_h5mu_path):
    output = random_h5mu_path()
    run_component(
        [
            "--input",
            input,
            "--output",
            output,
            "--output_compression",
            "gzip",
        ]
    )

    assert output.is_file(), "No output was created."
    data = mu.read_h5mu(output)

    assert "rna" in data.mod, 'Output should contain data.mod["rna"].'
    assert "X_tsne" in data.mod["rna"].obsm, "X_tsne should be present in .obsm"
    assert data.mod["rna"].obsm["X_tsne"].shape == (
        data.mod["rna"].n_obs,
        2,
    ), "t-SNE embedding should have shape (n_obs, 2)"


def test_obsm_output(run_component, random_h5mu_path):
    """Setting --obsm_output should write the embedding under that key."""
    output = random_h5mu_path()
    run_component(
        [
            "--input",
            input,
            "--output",
            output,
            "--output_compression",
            "gzip",
            "--obsm_output",
            "X_foo",
        ]
    )

    data = mu.read_h5mu(output)
    assert "X_foo" in data.mod["rna"].obsm, (
        "Embedding should be stored under the requested obsm key"
    )


def test_raise_if_obsm_input_missing(run_component, random_h5mu_path):
    """Should raise if the representation slot does not exist."""
    output = random_h5mu_path()
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input,
                "--output",
                output,
                "--output_compression",
                "gzip",
                "--obsm_input",
                "does_not_exist",
            ]
        )
    assert not output.is_file(), "No output should be created."
    assert (
        "ValueError: 'does_not_exist' was not found in .mod['rna'].obsm."
        in err.value.stdout.decode("utf-8")
    )


def test_overwrite_existing_slot(run_component, random_h5mu_path):
    """Writing into an existing .obsm slot should fail without --overwrite
    and succeed with it."""
    first = random_h5mu_path()
    run_component(["--input", input, "--output", first, "--output_compression", "gzip"])

    second = random_h5mu_path()
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            ["--input", str(first), "--output", second, "--output_compression", "gzip"]
        )
    assert "but field already exists" in err.value.stdout.decode("utf-8")

    run_component(
        [
            "--input",
            str(first),
            "--output",
            second,
            "--output_compression",
            "gzip",
            "--overwrite",
            "true",
        ]
    )
    assert "X_tsne" in mu.read_h5mu(second).mod["rna"].obsm


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
