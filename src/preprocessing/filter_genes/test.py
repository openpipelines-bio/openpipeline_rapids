import sys
import re
import pytest
import subprocess
import mudata as mu

## VIASH START
meta = {
    "name": "filter_genes",
    "resources_dir": "resources_test/",
    "config": "src/preprocessing/filter_genes/config.vsh.yaml",
    "executable": "target/docker/preprocessing/filter_genes/filter_genes",
}
## VIASH END

input = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"


def test_min_cells(run_component, random_h5mu_path):
    output = random_h5mu_path()
    run_component(
        [
            "--input",
            input,
            "--output",
            output,
            "--output_compression",
            "gzip",
            "--min_cells",
            "3",
        ]
    )

    assert output.is_file(), "No output was created."

    rna_in = mu.read_h5mu(input).mod["rna"]
    rna_out = mu.read_h5mu(output).mod["rna"]

    assert rna_out.n_obs == rna_in.n_obs, "Number of cells should be unchanged."
    assert rna_out.n_vars <= rna_in.n_vars, "Number of genes should not increase."
    assert rna_out.n_vars < rna_in.n_vars, "Some genes should have been filtered out."


def test_min_counts(run_component, random_h5mu_path):
    output = random_h5mu_path()
    run_component(
        [
            "--input",
            input,
            "--output",
            output,
            "--output_compression",
            "gzip",
            "--min_counts",
            "10",
        ]
    )

    assert output.is_file(), "No output was created."

    rna_in = mu.read_h5mu(input).mod["rna"]
    rna_out = mu.read_h5mu(output).mod["rna"]

    assert rna_out.n_obs == rna_in.n_obs, "Number of cells should be unchanged."
    assert rna_out.n_vars < rna_in.n_vars, "Some genes should have been filtered out."
    assert "n_counts" in rna_out.var.columns, (
        "filter_genes should annotate .var with n_counts."
    )
    assert "n_cells" in rna_out.var.columns, (
        "filter_genes should annotate .var with n_cells."
    )


def test_raise_if_no_threshold(run_component, random_h5mu_path):
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
            ]
        )
    assert re.search(
        r"Exactly one of --min_counts, --min_cells, --max_counts, --max_cells "
        r"must be set",
        err.value.stdout.decode("utf-8"),
    )


def test_raise_if_multiple_thresholds(run_component, random_h5mu_path):
    """rsc.pp.filter_genes accepts only one threshold per call."""
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
                "--min_counts",
                "10",
                "--max_counts",
                "10000",
            ]
        )
    assert re.search(
        r"Exactly one of --min_counts, --min_cells, --max_counts, --max_cells "
        r"must be set",
        err.value.stdout.decode("utf-8"),
    )


def test_layer(run_component, random_h5mu_path):
    """Filtering should run on the requested layer and keep it in the output."""
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
            "--min_cells",
            "3",
            "--layer",
            "counts",
        ]
    )

    rna_in = mu.read_h5mu(input).mod["rna"]
    rna_out = mu.read_h5mu(output).mod["rna"]
    assert rna_out.n_obs == rna_in.n_obs, "Number of cells should be unchanged."
    assert rna_out.n_vars < rna_in.n_vars, "Some genes should have been filtered out."
    assert "counts" in rna_out.layers, "The input layer should be preserved."


def test_raise_if_layer_missing(run_component, random_h5mu_path):
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
                "--min_cells",
                "3",
                "--layer",
                "does_not_exist",
            ]
        )
    assert re.search(
        r"Layer does_not_exist not found in modality rna",
        err.value.stdout.decode("utf-8"),
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
