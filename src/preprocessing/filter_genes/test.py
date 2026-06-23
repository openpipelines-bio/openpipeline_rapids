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
        r"At least one of --min_counts, --min_cells, --max_counts, --max_cells "
        r"must be set",
        err.value.stdout.decode("utf-8"),
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
