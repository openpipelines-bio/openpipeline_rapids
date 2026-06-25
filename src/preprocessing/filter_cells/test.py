import sys
import re
import pytest
import subprocess
import mudata as mu

## VIASH START
meta = {
    "name": "filter_cells",
    "resources_dir": "resources_test/",
    "config": "src/preprocessing/filter_cells/config.vsh.yaml",
    "executable": "target/docker/preprocessing/filter_cells/filter_cells",
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
        "--min_genes",
        "50",
    ]
    run_component(cmd_pars)

    assert output.is_file(), "No output was created."

    mu_input = mu.read_h5mu(input)
    mu_output = mu.read_h5mu(output)

    assert "rna" in mu_output.mod, 'Output should contain data.mod["rna"].'

    rna_in = mu_input.mod["rna"]
    rna_out = mu_output.mod["rna"]

    assert rna_out.n_vars == rna_in.n_vars, "Number of genes should be unchanged"
    assert rna_out.n_obs <= rna_in.n_obs, "Number of cells should not increase"
    assert rna_out.n_obs < rna_in.n_obs, "Some cells should have been filtered out"


def test_min_counts(run_component, random_h5mu_path):
    output = random_h5mu_path()
    cmd_pars = [
        "--input",
        input,
        "--output",
        output,
        "--output_compression",
        "gzip",
        "--min_counts",
        "1000",
    ]
    run_component(cmd_pars)

    assert output.is_file(), "No output was created."

    rna_in = mu.read_h5mu(input).mod["rna"]
    rna_out = mu.read_h5mu(output).mod["rna"]

    assert rna_out.n_vars == rna_in.n_vars, "Number of genes should be unchanged"
    assert rna_out.n_obs < rna_in.n_obs, "Some cells should have been filtered out"
    assert "n_counts" in rna_out.obs.columns, (
        "filter_cells should annotate .obs with n_counts."
    )
    assert "n_genes" in rna_out.obs.columns, (
        "filter_cells should annotate .obs with n_genes."
    )


def test_raise_if_no_threshold(run_component, random_h5mu_path):
    output = random_h5mu_path()
    cmd_pars = [
        "--input",
        input,
        "--output",
        output,
        "--output_compression",
        "gzip",
    ]
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(cmd_pars)
    assert re.search(
        r"Exactly one of --min_counts, --min_genes, --max_counts, --max_genes "
        r"must be set",
        err.value.stdout.decode("utf-8"),
    )


def test_raise_if_multiple_thresholds(run_component, random_h5mu_path):
    """rsc.pp.filter_cells accepts only one threshold per call."""
    output = random_h5mu_path()
    cmd_pars = [
        "--input",
        input,
        "--output",
        output,
        "--output_compression",
        "gzip",
        "--min_genes",
        "50",
        "--max_genes",
        "5000",
    ]
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(cmd_pars)
    assert re.search(
        r"Exactly one of --min_counts, --min_genes, --max_counts, --max_genes "
        r"must be set",
        err.value.stdout.decode("utf-8"),
    )


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
