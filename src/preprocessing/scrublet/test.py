import sys
import pytest
import subprocess
import mudata as mu

## VIASH START
meta = {
    "name": "scrublet",
    "resources_dir": "resources_test/",
    "config": "src/preprocessing/scrublet/config.vsh.yaml",
    "executable": "target/docker/preprocessing/scrublet/scrublet",
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

    assert rna_in.n_obs == rna_out.n_obs, (
        "Scrublet annotates cells, it should not filter them."
    )
    assert "doublet_score" in rna_out.obs, (
        'Output should contain .obs["doublet_score"].'
    )
    assert "predicted_doublet" in rna_out.obs, (
        'Output should contain .obs["predicted_doublet"].'
    )
    assert rna_out.obs["predicted_doublet"].dtype == bool, (
        "predicted_doublet should be a boolean column."
    )
    assert "scrublet" in rna_out.uns, 'Output should contain .uns["scrublet"].'


def test_random_state(run_component, random_h5mu_path):
    """Setting --random_state should still produce a valid annotated output."""
    output = random_h5mu_path()
    run_component(
        [
            "--input",
            input,
            "--output",
            output,
            "--output_compression",
            "gzip",
            "--random_state",
            "42",
        ]
    )

    assert output.is_file(), "No output was created."

    rna_out = mu.read_h5mu(output).mod["rna"]
    assert "doublet_score" in rna_out.obs, (
        'Output should contain .obs["doublet_score"].'
    )
    assert "predicted_doublet" in rna_out.obs, (
        'Output should contain .obs["predicted_doublet"].'
    )


def test_layer(run_component, random_h5mu_path):
    """Scrublet should run on the requested layer and leave .X untouched."""
    mu_in = mu.read_h5mu(input)
    mu_in.mod["rna"].layers["counts"] = mu_in.mod["rna"].X.copy()
    input_with_layer = random_h5mu_path()
    mu_in.write(input_with_layer)

    output = random_h5mu_path()
    run_component(
        [
            "--input",
            str(input_with_layer),
            "--output",
            output,
            "--output_compression",
            "gzip",
            "--layer",
            "counts",
        ]
    )

    assert output.is_file(), "No output was created."
    rna_out = mu.read_h5mu(output).mod["rna"]
    assert "doublet_score" in rna_out.obs, (
        'Output should contain .obs["doublet_score"].'
    )
    assert "predicted_doublet" in rna_out.obs, (
        'Output should contain .obs["predicted_doublet"].'
    )


def test_raise_if_layer_missing(run_component, random_h5mu_path):
    """Should raise if the requested layer does not exist."""
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
                "--layer",
                "does_not_exist",
            ]
        )
    assert not output.is_file(), "No output should be created."
    assert "Layer does_not_exist not found in modality rna" in err.value.stdout.decode(
        "utf-8"
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
