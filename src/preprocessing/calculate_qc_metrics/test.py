import sys
import re
import pytest
import subprocess
import mudata as mu

## VIASH START
meta = {
    "name": "calculate_qc_metrics",
    "resources_dir": "resources_test/",
    "config": "src/preprocessing/calculate_qc_metrics/config.vsh.yaml",
    "executable": "target/docker/preprocessing/calculate_qc_metrics/calculate_qc_metrics",
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

    data = mu.read_h5mu(output)
    assert "rna" in data.mod, 'Output should contain data.mod["rna"].'

    rna = data.mod["rna"]
    assert "n_genes_by_counts" in rna.obs.columns, (
        "Output .obs should contain n_genes_by_counts."
    )
    assert "total_counts" in rna.obs.columns, "Output .obs should contain total_counts."
    assert "n_cells_by_counts" in rna.var.columns, (
        "Output .var should contain n_cells_by_counts."
    )
    assert "total_counts" in rna.var.columns, "Output .var should contain total_counts."


def test_qc_vars(run_component, random_h5mu_path):
    """Providing a boolean .var column via --qc_vars should add per-cell
    percentage metrics for that group of genes."""
    mu_orig = mu.read_h5mu(input)
    rna = mu_orig.mod["rna"]
    rna.var["mito"] = rna.var_names.str.startswith("MT-")
    input_with_qc_var = random_h5mu_path()
    mu_orig.write(input_with_qc_var)

    output = random_h5mu_path()
    run_component(
        [
            "--input",
            str(input_with_qc_var),
            "--output",
            output,
            "--output_compression",
            "gzip",
            "--qc_vars",
            "mito",
        ]
    )

    assert output.is_file(), "No output was created."

    rna_out = mu.read_h5mu(output).mod["rna"]
    assert "pct_counts_mito" in rna_out.obs.columns, (
        "Output .obs should contain pct_counts_mito."
    )
    assert "total_counts_mito" in rna_out.obs.columns, (
        "Output .obs should contain total_counts_mito."
    )


def test_log1p(run_component, random_h5mu_path):
    """--log1p adds log1p-transformed metrics on top of the regular ones; with
    --log1p false only the regular metrics are present."""
    with_log1p = random_h5mu_path()
    run_component(
        [
            "--input",
            input,
            "--output",
            with_log1p,
            "--output_compression",
            "gzip",
            "--log1p",
            "true",
        ]
    )
    rna_with = mu.read_h5mu(with_log1p).mod["rna"]
    assert "total_counts" in rna_with.obs.columns, (
        "Regular metrics should always be present."
    )
    assert "log1p_total_counts" in rna_with.obs.columns, (
        "--log1p should add log1p-transformed .obs metrics on top of the regular ones."
    )
    assert "log1p_total_counts" in rna_with.var.columns, (
        "--log1p should add log1p-transformed .var metrics on top of the regular ones."
    )

    without_log1p = random_h5mu_path()
    run_component(
        [
            "--input",
            input,
            "--output",
            without_log1p,
            "--output_compression",
            "gzip",
            "--log1p",
            "false",
        ]
    )
    rna_without = mu.read_h5mu(without_log1p).mod["rna"]
    assert "total_counts" in rna_without.obs.columns, (
        "Regular metrics should still be present without --log1p."
    )
    assert "log1p_total_counts" not in rna_without.obs.columns, (
        "No log1p metrics should be added when --log1p is false."
    )
    assert "log1p_total_counts" not in rna_without.var.columns, (
        "No log1p metrics should be added when --log1p is false."
    )


def test_raise_if_layer_missing(run_component, random_h5mu_path):
    """Requesting a layer that does not exist should fail."""
    output = random_h5mu_path()
    with pytest.raises(subprocess.CalledProcessError) as err:
        run_component(
            [
                "--input",
                input,
                "--output",
                output,
                "--layer",
                "nonexistent",
            ]
        )
    assert re.search(
        r"Layer 'nonexistent' not found in modality 'rna'",
        err.value.stdout.decode("utf-8"),
    )


def test_empty_layer_uses_x(run_component, random_h5mu_path):
    """An empty --layer means ".X"; it must not be treated as a layer name."""
    output = random_h5mu_path()
    run_component(
        [
            "--input",
            input,
            "--output",
            output,
            "--layer",
            "",
        ]
    )

    assert output.is_file(), "No output was created."


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
