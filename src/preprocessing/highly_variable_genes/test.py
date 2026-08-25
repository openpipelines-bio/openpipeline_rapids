import sys
import pytest
import mudata as mu
import numpy as np
import scanpy as sc
from openpipeline_testutils.asserters import assert_annotation_objects_equal

## VIASH START
meta = {
    "name": "highly_variable_genes",
    "resources_dir": "resources_test/",
    "config": "src/preprocessing/highly_variable_genes/config.vsh.yaml",
    "executable": "target/docker/preprocessing/highly_variable_genes/highly_variable_genes",
}
## VIASH END

input = f"{meta['resources_dir']}/pbmc_1k_protein_v3/pbmc_1k_protein_v3_filtered_feature_bc_matrix.h5mu"


@pytest.fixture
def lognormed_input(random_h5mu_path):
    """rapids-singlecell hvg with the default flavor expects log-normalized
    data. sc.pp.log1p preserves sparsity."""
    mu_in = mu.read_h5mu(input)
    rna = mu_in.mod["rna"]
    sc.pp.log1p(rna)
    rna.layers["log_normalized"] = rna.X.copy()
    path = random_h5mu_path()
    mu_in.write(path)
    return path


@pytest.fixture
def lognormed_filtered_input(random_h5mu_path):
    """cell_ranger flavor requires gene-filtered data - running it on the
    raw filtered_feature_bc_matrix collapses pandas.cut's bin edges to
    identical near-zero means."""
    mu_in = mu.read_h5mu(input)
    rna = mu_in.mod["rna"]
    sc.pp.filter_genes(rna, min_counts=20)
    sc.pp.log1p(rna)
    rna.layers["log_normalized"] = rna.X.copy()
    path = random_h5mu_path()
    mu_in.write(path)
    return path


def test_run(run_component, random_h5mu_path, lognormed_input):
    output = random_h5mu_path()
    cmd_pars = [
        "--input",
        lognormed_input,
        "--output",
        output,
        "--output_compression",
        "gzip",
    ]
    run_component(cmd_pars)

    assert output.is_file(), "No output was created."

    mu_input = mu.read_h5mu(lognormed_input)
    mu_output = mu.read_h5mu(output)

    assert "rna" in mu_output.mod, 'Output should contain data.mod["rna"].'

    rna_in = mu_input.mod["rna"]
    rna_out = mu_output.mod["rna"]

    assert rna_in.shape == rna_out.shape, "Should have same shape as before"

    # Default --var_name_filter is "highly_variable"
    assert "highly_variable" in rna_out.var.columns, (
        "Output should contain a 'highly_variable' boolean column in .var"
    )
    assert rna_out.var["highly_variable"].dtype == bool
    assert rna_out.var["highly_variable"].any(), "At least one HVG should be selected"

    # Default --varm_name is "hvg"
    assert "hvg" in rna_out.varm, "Output should contain an 'hvg' .varm slot"
    assert rna_out.varm["hvg"].shape[0] == rna_out.n_vars
    # The per-gene dispersion metrics should land in .varm (openpipeline#1143)
    assert {
        "highly_variable",
        "means",
        "dispersions",
        "dispersions_norm",
    }.issubset(rna_out.varm["hvg"].columns), (
        "Dispersion metrics should be stored in the 'hvg' .varm slot"
    )
    assert rna_out.varm["hvg"]["highly_variable"].equals(
        rna_out.var["highly_variable"]
    ), ".varm and .var should agree on which features are highly variable"

    # rapids-singlecell records the flavor in .uns (openpipeline#1141)
    assert rna_out.uns["hvg"]["flavor"] == "seurat", (
        "The flavor should be recorded in .uns['hvg']"
    )

    # X should be untouched - hvg only annotates .var / .varm / .uns
    assert np.allclose(rna_out.X.toarray(), rna_in.X.toarray()), (
        "X should be untouched by highly_variable_genes"
    )

    # Strip the new annotations and compare the rest of the object
    rna_out.var = rna_out.var.drop(columns=["highly_variable"])
    del rna_out.varm["hvg"]
    del rna_out.uns["hvg"]
    if "rna:highly_variable" in mu_output.var.columns:
        mu_output.var = mu_output.var.drop(columns=["rna:highly_variable"])
    assert_annotation_objects_equal(mu_input, mu_output)


def test_n_top_features(run_component, random_h5mu_path, lognormed_input):
    """Setting --n_top_features should produce exactly that many HVGs."""
    output = random_h5mu_path()
    run_component(
        [
            "--input",
            lognormed_input,
            "--output",
            output,
            "--output_compression",
            "gzip",
            "--n_top_features",
            "50",
        ]
    )

    rna_out = mu.read_h5mu(output).mod["rna"]
    assert int(rna_out.var["highly_variable"].sum()) == 50, (
        "n_top_features should yield exactly that many highly variable genes"
    )


def test_flavor_cell_ranger(run_component, random_h5mu_path, lognormed_filtered_input):
    """Cell Ranger flavor should give a different HVG set than Seurat."""
    output_seurat = random_h5mu_path()
    output_cr = random_h5mu_path()

    run_component(
        [
            "--input",
            lognormed_filtered_input,
            "--output",
            output_seurat,
            "--output_compression",
            "gzip",
            "--n_top_features",
            "100",
            "--flavor",
            "seurat",
        ]
    )
    run_component(
        [
            "--input",
            lognormed_filtered_input,
            "--output",
            output_cr,
            "--output_compression",
            "gzip",
            "--n_top_features",
            "100",
            "--flavor",
            "cell_ranger",
        ]
    )

    seurat_hvg = mu.read_h5mu(output_seurat).mod["rna"].var["highly_variable"]
    cr_hvg = mu.read_h5mu(output_cr).mod["rna"].var["highly_variable"]

    assert int(seurat_hvg.sum()) == 100
    assert int(cr_hvg.sum()) == 100
    assert not seurat_hvg.equals(cr_hvg), (
        "Different flavors should select different highly variable genes"
    )


def test_flavor_seurat_v3_requires_n_top_features(run_component, random_h5mu_path):
    """flavor='seurat_v3' without --n_top_features must error out."""
    import subprocess

    output = random_h5mu_path()
    with pytest.raises(subprocess.CalledProcessError):
        run_component(
            [
                "--input",
                input,
                "--output",
                output,
                "--flavor",
                "seurat_v3",
            ]
        )


def test_flavor_seurat_v3(run_component, random_h5mu_path):
    """flavor='seurat_v3' expects raw counts and requires --n_top_features."""
    output = random_h5mu_path()
    run_component(
        [
            "--input",
            input,
            "--output",
            output,
            "--output_compression",
            "gzip",
            "--flavor",
            "seurat_v3",
            "--n_top_features",
            "50",
        ]
    )

    rna_out = mu.read_h5mu(output).mod["rna"]
    assert int(rna_out.var["highly_variable"].sum()) == 50
    # seurat_v3 stores variance-based metrics; ensure they reach .varm
    # rather than being dropped (openpipeline#1143)
    assert {
        "highly_variable",
        "means",
        "variances",
        "variances_norm",
        "highly_variable_rank",
    }.issubset(rna_out.varm["hvg"].columns), (
        "seurat_v3 variance metrics should be stored in the 'hvg' .varm slot"
    )
    assert int(rna_out.varm["hvg"]["highly_variable"].sum()) == 50
    assert rna_out.uns["hvg"]["flavor"] == "seurat_v3", (
        "The flavor should be recorded in .uns['hvg']"
    )


def test_input_layer(run_component, random_h5mu_path):
    """Using --input_layer should give the same result as running on a mudata
    whose X is that same matrix."""
    mu_in = mu.read_h5mu(input)
    rna = mu_in.mod["rna"]
    sc.pp.log1p(rna)
    rna.layers["log_normalized"] = rna.X.copy()
    layer_input = random_h5mu_path()
    mu_in.write(layer_input)

    output = random_h5mu_path()
    run_component(
        [
            "--input",
            str(layer_input),
            "--output",
            output,
            "--output_compression",
            "gzip",
            "--input_layer",
            "log_normalized",
            "--n_top_features",
            "100",
        ]
    )

    rna_out = mu.read_h5mu(output).mod["rna"]
    assert "highly_variable" in rna_out.var.columns
    assert int(rna_out.var["highly_variable"].sum()) == 100
    # The layer should still be present
    assert "log_normalized" in rna_out.layers


def test_var_name_filter_and_varm_name(
    run_component, random_h5mu_path, lognormed_input
):
    """Custom --var_name_filter and --varm_name should rename the output slots."""
    output = random_h5mu_path()
    run_component(
        [
            "--input",
            lognormed_input,
            "--output",
            output,
            "--output_compression",
            "gzip",
            "--var_name_filter",
            "filter_with_hvg",
            "--varm_name",
            "hvg_metrics",
            "--n_top_features",
            "50",
        ]
    )

    rna_out = mu.read_h5mu(output).mod["rna"]
    assert "filter_with_hvg" in rna_out.var.columns
    assert "highly_variable" not in rna_out.var.columns
    assert "hvg_metrics" in rna_out.varm
    assert "hvg" not in rna_out.varm
    assert int(rna_out.var["filter_with_hvg"].sum()) == 50
    # .varm keeps the flag under the name rapids-singlecell assigned it
    assert "highly_variable" in rna_out.varm["hvg_metrics"].columns
    assert rna_out.varm["hvg_metrics"]["highly_variable"].equals(
        rna_out.var["filter_with_hvg"].rename("highly_variable")
    )


def test_obs_batch_key(run_component, random_h5mu_path, lognormed_input):
    """With --obs_batch_key, the per-batch HVG metadata columns should appear."""
    mu_in = mu.read_h5mu(lognormed_input)
    rna = mu_in.mod["rna"]
    rna.obs["batch"] = "A"
    column_index = rna.obs.columns.get_indexer(["batch"])
    rna.obs.iloc[slice(rna.n_obs // 2, None), column_index] = "B"
    batch_input = random_h5mu_path()
    mu_in.write(batch_input)

    output = random_h5mu_path()
    run_component(
        [
            "--input",
            str(batch_input),
            "--output",
            output,
            "--output_compression",
            "gzip",
            "--obs_batch_key",
            "batch",
            "--n_top_features",
            "100",
        ]
    )

    rna_out = mu.read_h5mu(output).mod["rna"]
    assert "highly_variable" in rna_out.var.columns
    # batched mode emits these per-batch columns into .varm["hvg"]
    hvg_metrics = rna_out.varm["hvg"]
    assert "highly_variable_nbatches" in hvg_metrics.columns
    assert "highly_variable_intersection" in hvg_metrics.columns


def test_empty_input_layer_uses_x(run_component, random_h5mu_path, lognormed_input):
    """An empty --input_layer means ".X"; it must not be treated as a layer name."""
    output = random_h5mu_path()
    run_component(
        [
            "--input",
            lognormed_input,
            "--output",
            output,
            "--input_layer",
            "",
        ]
    )

    assert output.is_file(), "No output was created."


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
