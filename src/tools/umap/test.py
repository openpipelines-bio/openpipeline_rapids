import sys
import pytest
import subprocess
import mudata as mu

## VIASH START
meta = {
    "name": "umap",
    "resources_dir": "resources_test/",
    "config": "src/tools/umap/config.vsh.yaml",
    "executable": "target/docker/tools/umap/umap",
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
    assert "X_umap" in data.mod["rna"].obsm, "X_umap should be present in .obsm"
    assert data.mod["rna"].obsm["X_umap"].shape == (
        data.mod["rna"].n_obs,
        2,
    ), "UMAP embedding should have shape (n_obs, 2)"


def test_num_components(run_component, random_h5mu_path):
    """Setting --num_components should change the embedding dimensionality."""
    output = random_h5mu_path()
    run_component(
        [
            "--input",
            input,
            "--output",
            output,
            "--output_compression",
            "gzip",
            "--num_components",
            "5",
        ]
    )

    data = mu.read_h5mu(output)
    assert data.mod["rna"].obsm["X_umap"].shape == (
        data.mod["rna"].n_obs,
        5,
    ), "UMAP embedding should have the requested number of components"


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


def test_raise_if_uns_neighbors_missing(run_component, random_h5mu_path):
    """Should raise if the neighbors slot does not exist."""
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
                "--uns_neighbors",
                "does_not_exist",
            ]
        )
    assert not output.is_file(), "No output should be created."
    assert (
        "ValueError: 'does_not_exist' was not found in .mod['rna'].uns."
        in err.value.stdout.decode("utf-8")
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
