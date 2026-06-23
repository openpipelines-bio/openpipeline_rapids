import sys
import pytest
import subprocess
import mudata as mu

## VIASH START
meta = {
    "name": "harmony_integrate",
    "resources_dir": "resources_test/",
    "config": "src/preprocessing/harmony_integrate/config.vsh.yaml",
    "executable": "target/docker/preprocessing/harmony_integrate/harmony_integrate",
}
## VIASH END

input = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_mms.h5mu"


def _prepare_input_with_batch(random_h5mu_path):
    """Read the input (which already has .obsm['X_pca']), add a categorical
    batch column to .obs, and write it back to a fresh h5mu path."""
    mdata = mu.read_h5mu(input)
    rna = mdata.mod["rna"]
    rna.obs["batch"] = ["a" if i % 2 else "b" for i in range(rna.n_obs)]
    path = random_h5mu_path()
    mdata.write(path)
    return path


def test_run(run_component, random_h5mu_path):
    input_with_batch = _prepare_input_with_batch(random_h5mu_path)
    output = random_h5mu_path()
    run_component(
        [
            "--input",
            str(input_with_batch),
            "--output",
            output,
            "--output_compression",
            "gzip",
            "--obs_covariates",
            "batch",
        ]
    )

    assert output.is_file(), "No output was created."

    data = mu.read_h5mu(output)
    assert "rna" in data.mod, 'Output should contain data.mod["rna"].'

    rna_out = data.mod["rna"]
    assert "X_pca_harmony" in rna_out.obsm, "Output should have .obsm['X_pca_harmony']"
    assert rna_out.obsm["X_pca_harmony"].shape == rna_out.obsm["X_pca"].shape, (
        "Integrated embedding should have the same shape as the input PCA."
    )


def test_custom_obsm_output(run_component, random_h5mu_path):
    """A custom --obsm_output slot is stored under exactly the requested key."""
    input_with_batch = _prepare_input_with_batch(random_h5mu_path)
    output = random_h5mu_path()
    run_component(
        [
            "--input",
            str(input_with_batch),
            "--output",
            output,
            "--output_compression",
            "gzip",
            "--obs_covariates",
            "batch",
            "--obsm_output",
            "X_custom_harmony",
        ]
    )

    rna_out = mu.read_h5mu(output).mod["rna"]
    assert "X_custom_harmony" in rna_out.obsm, (
        "Output should have .obsm['X_custom_harmony']"
    )
    assert "X_pca_harmony" not in rna_out.obsm


def test_raise_if_obsm_input_missing(run_component, random_h5mu_path):
    """A non-existent --obsm_input should raise an error."""
    input_with_batch = _prepare_input_with_batch(random_h5mu_path)
    output = random_h5mu_path()
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                str(input_with_batch),
                "--output",
                output,
                "--output_compression",
                "gzip",
                "--obs_covariates",
                "batch",
                "--obsm_input",
                "X_does_not_exist",
            ]
        )
    assert "obsm slot 'X_does_not_exist' not found" in err.value.stdout.decode("utf-8")


def test_raise_if_covariate_missing(run_component, random_h5mu_path):
    """A non-existent --obs_covariates column should raise an error."""
    input_with_batch = _prepare_input_with_batch(random_h5mu_path)
    output = random_h5mu_path()
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                str(input_with_batch),
                "--output",
                output,
                "--output_compression",
                "gzip",
                "--obs_covariates",
                "not_a_column",
            ]
        )
    assert "obs column 'not_a_column' not found" in err.value.stdout.decode("utf-8")


def test_overwrite_existing_slot(run_component, random_h5mu_path):
    """Writing into an existing .obsm slot should fail without --overwrite
    and succeed with it."""
    input_with_batch = _prepare_input_with_batch(random_h5mu_path)
    first = random_h5mu_path()
    run_component(
        [
            "--input",
            str(input_with_batch),
            "--output",
            first,
            "--output_compression",
            "gzip",
            "--obs_covariates",
            "batch",
        ]
    )

    second = random_h5mu_path()
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                str(first),
                "--output",
                second,
                "--output_compression",
                "gzip",
                "--obs_covariates",
                "batch",
            ]
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
            "--obs_covariates",
            "batch",
            "--overwrite",
            "true",
        ]
    )
    assert "X_pca_harmony" in mu.read_h5mu(second).mod["rna"].obsm


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
